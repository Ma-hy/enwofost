# coding: utf-8

import os,sys
import pcse
import datetime as dt
import numpy as np
import gdal
from netCDF4 import Dataset,date2index
import os.path
import cdsapi
from pcse.fileinput import CABOFileReader
from pcse.base_classes import ParameterProvider
from pcse.util import WOFOST71SiteDataProvider
from pcse.fileinput import YAMLAgroManagementReader
from pcse.fileinput import CABOWeatherDataProvider
from pcse.models import Wofost71_PP
import scipy.stats as ss
import copy
from bisect import bisect
from pcse.db import NASAPowerWeatherDataProvider

home = os.path.dirname(os.path.realpath("__file__"))
data_dir = home+"/data/"

def calculate_hum(tdew):
    from numpy import exp
    tdew=tdew-273.15
    tmp = (17.27 * tdew) / (tdew + 237.3)
    ea = 0.6108 * exp(tmp)
    return ea

def retrieve_pixel_value(geo_coord, data_source):
    """Return floating-point value that corresponds to given point."""
    dataset = gdal.Open(data_source)
    if dataset is None:
        raise Exception("Wrong data source! Please check it.")
    band = dataset.GetRasterBand(1)

    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    transform = dataset.GetGeoTransform()

    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    data = band.ReadAsArray(0, 0, cols, rows)
    col = int((geo_coord[0] - xOrigin) / pixelWidth)
    row = int((yOrigin - geo_coord[1] ) / pixelHeight)
    dataset = None

    return data[row][col]
    
def gen_era_cabo(mylat, mylon, start_year, end_year, inputfile=data_dir,
                 data_dir=None):
    size = 0.25
    station_number=1
    site="%5.2f_%5.2f"%(int((mylon+size/2.)/size)*size,int((mylat+size/2.)/size)*size)
    c1=-0.18
    c2=-0.55

    elev=retrieve_pixel_value([mylon,mylat],data_dir+"alwdgg.tif")
    parnames = ["ssrd","mx2t","mn2t","tp","u10","v10","d2m"]
    createVar = globals()

    for year in range(start_year,end_year+1):
        fname = inputfile+site+".%s"%(str(year)[-3:])
        if not os.path.isfile(fname):  # if CABO weather file doesn't exist
            fname = inputfile + "%d_%d_10d_%s.nc"%(int(mylon/10.)*10,int(mylat/10.)*10,str(year))
            #fname = inputfile+"%s_%s_%s.nc"%(int(mylon),int(mylat),str(year))
            if not os.path.isfile(fname):  # if NetCDF weather file doesn't exist
                print("There's no weather driver avaliable, you need to download it.")
                query_nc = input("Do you want to download from ECMWF? It may take more than one hour: (y/n)")
                if query_nc == 'y'or 'Y':
                    c = cdsapi.Client()                
                    c.retrieve(
                        'reanalysis-era5-single-levels',
                        {
                            'variable':[
                                'surface_solar_radiation_downwards',
                                'maximum_2m_temperature_since_previous_post_processing','minimum_2m_temperature_since_previous_post_processing',
                                'total_precipitation',
                                '10m_u_component_of_wind','10m_v_component_of_wind',
                                '2m_dewpoint_temperature'                            
                            ],
                            'product_type':'reanalysis',
                            'year':year,
                            'month':[
                                '01','02','03',
                                '04','05','06',
                                '07','08','09',
                                '10','11','12'
                            ],
                            'day':[
                                '01','02','03',
                                '04','05','06',
                                '07','08','09',
                                '10','11','12',
                                '13','14','15',
                                '16','17','18',
                                '19','20','21',
                                '22','23','24',
                                '25','26','27',
                                '28','29','30',
                                '31'
                            ],
                            'time':[
                                '00:00','01:00','02:00',
                                '03:00','04:00','05:00',
                                '06:00','07:00','08:00',
                                '09:00','10:00','11:00',
                                '12:00','13:00','14:00',
                                '15:00','16:00','17:00',
                                '18:00','19:00','20:00',
                                '21:00','22:00','23:00'
                            ],
                            'area' : "%d/%d/%d/%d"%(int(mylat/10.+1.)*10,int(mylon/10.)*10,int(mylat/10.)*10,int(mylon/10.+1.)*10,),
                            'format':'netcdf'
                        },
                        fname)
            print("Generating cabo format weather driver from a ERA5 NetCDF file: ",fname)
            dataset=Dataset(fname)
            for par in parnames:    
                createVar[par]=dataset.variables[par][:]
            uplat, dnlat = dataset.variables["latitude"][:].max(),dataset.variables["latitude"][:].min()
            uplon, dnlon = dataset.variables["longitude"][:].max(),dataset.variables["longitude"][:].min()
            x=int((mylon-dnlon+size/2)/size)
            y=int((mylat-uplat-size/2)/-size)
            times=dataset.variables["time"]
                
            rad = np.sum(ssrd.reshape(-1, 24,ssrd.shape[1],ssrd.shape[2]), axis=1)/1000.
            tmax = np.max(mx2t.reshape(-1, 24,mx2t.shape[1],mx2t.shape[2]), axis=1)-273.15
            tmin = np.min(mn2t.reshape(-1, 24,mn2t.shape[1],mn2t.shape[2]), axis=1)-273.15
            prec = np.sum(tp.reshape(-1, 24,tp.shape[1],tp.shape[2]), axis=1)*1000.
            prec[prec < 0.01] = 0
            wind_u = np.mean(u10.reshape(-1,24,u10.shape[1],u10.shape[2]),axis=1)
            wind_v = np.mean(v10.reshape(-1,24,v10.shape[1],v10.shape[2]),axis=1)
            wind = np.sqrt(np.square(wind_u)+np.square(wind_v))
            hum = calculate_hum(np.mean(d2m.reshape(-1, 24,d2m.shape[1],d2m.shape[2]), axis=1))

            f=open(inputfile+site+".%s"%(str(year)[-3:]),"w+")
            print (year,"begins for site: ","%s in the year of %s"%(site,str(year)))
            f.write("*------------------------------------------------------------*"+"\n"
                        +'*'+"%12s"%("Country: ")+"By Coordinate"+"\n"
                        +'*'+"%12s"%("Station: ")+"%s"%site+"\n"
                        +'*'+"%12s"%("Year: ")+"%d"%(year)+"\n"
                        +'*'+"%12s"%("Origin: ")+"ERA5 Reanalysis"+"\n"
                        +'*'+"%12s"%("Author: ")+"UCL, MaHongYuan"+"\n"
                        +'*'+"%12s"%("Longitude: ")+"%f"%(mylon)+" E"+"\n"
                        +'*'+"%12s"%("Latitude: ")+"%f"%(mylat)+" N"+"\n"
                        +'*'+"%12s"%("Elevation: ")+"%.2f"%(elev)+" m"+"\n"
                        +'*'+"%12s"%("Columns: ")+"\n"
                        +'*'+"%12s"%("======== ")+"\n"
                        +'*'+"  station number"+"\n"
                        +'*'+"  year"+"\n"
                        +'*'+"  day"+"\n"
                        +'*'+"  irradiation (kJ·m-2·d-1)"+"\n"
                        +'*'+"  minimum temperature (degrees Celsius)"+"\n"
                        +'*'+"  maximum temperature (degrees Celsius)"+"\n"
                        +'*'+"  vapour pressure (kPa)"+"\n" 
                        +'*'+"  mean wind speed (m·s-1)"+"\n" 
                        +'*'+"  precipitation (mm·d-1)"+"\n" 
                        +'**'+" WCCDESCRIPTION="+site+", China"+"\n" 
                        +'**'+" WCCFORMAT=2"+"\n" 
                        +'**'+" WCCYEARNR="+"%d"%(year)+"\n" 
                        +"*------------------------------------------------------------*"+"\n"
                        +"%.2f  %.2f  %.2f  %.2f  %.2f\n"%(mylon, mylat, elev, c1, c2)
                        )

            year_s = date2index(dt.datetime.strptime('%s-01-01 00:00:00'%year, '%Y-%m-%d %H:%M:%S'),times)
            year_e = date2index(dt.datetime.strptime('%s-12-31 00:00:00'%year, '%Y-%m-%d %H:%M:%S'),times)
            for d in range(rad.shape[0]):
                f.write("%d"%(station_number)+"\t"+"%d"%(year)+"\t"+"%3d"%(1+d)+"\t"
                            +"%5d"%(round(rad[d,y,x]))+"\t"
                            +"%5.1f"%(round(tmin[d,y,x]*10)/10)+"\t"
                            +"%5.1f"%(round(tmax[d,y,x]*10)/10)+"\t"
                            +"%5.3f"%(round(hum[d,y,x]*1000)/1000)+"\t"
                            +"%4.1f"%(round(wind[d,y,x]*10)/10)+"\t"
                            +"%4.1f"%(round(prec[d,y,x]*10)/10)+"\n")
            f.close()    
            
    dataset = None    


def define_prior_distributions(chunk=data_dir+"par_prior.csv"):
    """A function to interpret the prior distributions from an Excel c&p job.
    Returns a dictionary indexed by parameter name and a pointer to a
    scipy.stats function that provides a logpdf method ;-). It also returns a
    list that can be used to map from vector to dictionary."""
    prior_dist = {}
    param_list = []
    param_xvalue = {}
    param_type = {}
    if not os.path.isfile(chunk):
        raise ValueError("There is no valid configuration file")
    with open(chunk,'r') as chunkfile:
        chunk_data = chunkfile.readlines()
    for line in chunk_data:
        if not line.strip().startswith("#"):
            if line.find("Uniform") >= 0:
                param,ty, xp, yp, xmin, xmax, sigma, _ = line.split(",")
                xmin = float(xmin)
                xmax = float(xmax)               
                dist = ss.uniform(xmin, (xmax - xmin))
            elif line.find("Gaussian") >= 0:
                param,ty, xp, yp, xmin, xmax, sigma, _ = line.split(",")
                if ty[0] == "X":
                    xp,yp = yp,xp
                lower = float(xmin)
                upper = float(xmax)
                mu = float(yp)
                sigma = float(sigma)
                dist = ss.truncnorm((lower - mu) / sigma,
                                    (upper - mu) / sigma, loc=mu, scale=sigma)            
            prior_dist[param] = dist
            param_list.append(param)
            param_xvalue[param] = float(xp)
            param_type[param] = ty
    return prior_dist, param_list, param_xvalue,param_type

def ensemble_wofost(lon = 115.55, lat=38., start = dt.date(2008,10,12),
                    end = None, en_size = 3, prior_file = None,
                    weather_type = "NASA", weather_path = None, out_en_file = None, data_dir=None):
    """
    This is a function to generate a emsemble of WOFOST paramters and corresponding output.
    you need to specify Longitude (lon), Latitude (lat), 
    start time of crop (start), end time of crop (end, 270 days duration by default),
    emsemble size (en_size), configuration file for prior distributions of pramaters (prior_file), 
    weather driver dataset type (weather_type), it's set to NASA Power dataset "NASA" by default,
    you could use ERA5 "ERA5" instead or use your own CABO file (%your_cabo_files_name%).)
    """
    if data_dir is None:
        #home = os.path.dirname(os.path.realpath("__file__"))
        os.path.split(os.path.realpath(__file__))[0]
        data_dir = home+"/data/"
    if prior_file is None:
        prior_file = data_dir+"par_prior.csv"    
    if out_en_file is None:
        out_en_file = data_dir+"WOFOST_par_ensemble.npy"
    
    if lat < -90 or lat > 90:
        msg = "Latitude should be between -90 and 90 degrees."
        raise ValueError(msg)
    if lon < -180 or lon > 180:
        msg = "Longitude should be between -180 and 180 degrees."
        raise ValueError(msg)
    if end == None:
        end = start + dt.timedelta(days=270)
    if start >= end:
        msg = "Start time should be earlier than end time."
        raise ValueError(msg)
    if weather_type == "NASA":
        print("Downloading weather driver from NASA Power...")
        weather = NASAPowerWeatherDataProvider(latitude=lat, longitude=lon)
    elif weather_type[:3] == "ERA" or weather_type[:3] == "era":
        print("ERA5 reanalysis dataset used.")
        if  weather_path is None or not os.path.isdir(weather_path):
            msg = "Please provide a valid path for weahter driver data."
            raise ValueError(msg)
        gen_era_cabo(lat, lon, start.year, end.year, inputfile=weather_path, 
                     data_dir=data_dir)
        size = 0.25
        weather_name = "%5.2f_%5.2f"%(int((lon+size/2.)/size)*size,int((lat+size/2.)/size)*size)
        weather = CABOWeatherDataProvider(weather_name, fpath=weather_path)
    else:
        if weather_path == None:
            raise ValueError("Please provide your weather driver path!")
        weather = CABOWeatherDataProvider(weather_type, fpath=weather_path)
        
        
    varnames = ["day", "TAGP", "LAI", "TWSO","DVS"]
    tmp={}
    
    cropfile = os.path.join(data_dir, 'WWH108.CAB')
    crop = CABOFileReader(cropfile)
    soilfile = os.path.join(data_dir, 'Hengshui.soil')
    soil = CABOFileReader(soilfile)
    site = WOFOST71SiteDataProvider(WAV=100, CO2=360)
    parameters = ParameterProvider(soildata=soil, cropdata=crop, sitedata=site)
    agromanagement_file = os.path.join(data_dir, 'shenzhou_wheat.amgt')
    agromanagement = YAMLAgroManagementReader(agromanagement_file)
    (key, value), = agromanagement[0].items()
    agromanagement[0][start]=agromanagement[0].pop(key)
    value['CropCalendar']['crop_start_date']=start
    value['CropCalendar']['crop_end_date']=end
       
    prior_dist, prior_list, param_xvalue,param_type = define_prior_distributions(chunk=prior_file)
    z_start = np.empty((len(prior_list), en_size))
    for i, param in enumerate(prior_list):
        z_start[i, :] = prior_dist[param].rvs(en_size)
    outdata = []

    for i in range(en_size):
        theta_dict = dict(zip(prior_list, z_start[:,i]))
        cropdata = copy.deepcopy(crop)
        tb_x = {}
        tb_y = {}
        tb_t = {}
        tmp_dict={}
        for par in theta_dict.keys():
            try:
                if param_type[par] != 'S':
                    tb_index = par.find("TB")
                    if tb_index < 0:
                        print(param_xvparam_typealue[par])
                        raise Exception("Are you sure %s is a table value?"%par)
                    tb_name = par[:tb_index+2]
                    tmp_list = [param_xvalue[par],theta_dict[par]]
                    if not tb_name in tb_x:
                        tb_x[tb_name] = np.array([param_xvalue[par]])
                        tb_y[tb_name] = np.array([theta_dict[par]])
                        tb_t[tb_name] = param_type[par]
                    else:
                        tb_x[tb_name] = np.append(tb_x[tb_name],param_xvalue[par])
                        tb_y[tb_name] = np.append(tb_y[tb_name],theta_dict[par])
            except KeyError:
                raise Exception("There's something wrong with %s, please check it."%par)
        tmp_dict={}
        for par in tb_x.keys():  # Table parameters
            s_i = np.argsort(tb_x[par])
            s_x = tb_x[par][s_i]
            s_v = tb_y[par][s_i]
            par_tb = []
    #         print(par,tb_t[par],cropdata[par],s_x,s_v)
            if tb_t[par][1] == 'P':   
                for i in range(len(tb_x[par])):
                    if tb_t[par][0] == 'Y':               # Partly change table Y values
                        if s_x[i] in cropdata[par][::2]:  # change old value
                            c_i = cropdata[par][::2].index(s_x[i])
                            cropdata[par][c_i*2] = s_v[i]
                        else:                             # insert new value
                            array_X = cropdata[par][::2]
                            array_Y = cropdata[par][1:][::2]
                            ins_i = bisect(array_X, s_x[i])
                            cropdata[par].insert( ins_i*2, s_x[i])
                            cropdata[par].insert( ins_i*2+1, s_v[i])
                        #print(cropdata[par])
                    else:                                 # Partly change table X values
                        if s_x[i] in cropdata[par][1:][::2]:  # change old value
                            c_i = cropdata[par][1:][::2].index(s_x[i])
                            cropdata[par][c_i*2] = s_v[i]
                        else:                             # insert new value
                            array_X = cropdata[par][::2]
                            array_Y = cropdata[par][1:][::2]
                            ins_i = bisect(array_X, s_x[i])
                            cropdata[par].insert( ins_i*2, s_x[i])
                            cropdata[par].insert( ins_i*2+1, s_v[i])
                        #print(cropdata[par])
            elif tb_t[par][1] == 'A':                     
                if tb_t[par][0] == 'Y':                  # Totally change table Y values
                    for i in range(len(tb_x[par])):
                        par_tb.append(s_x[i])
                        par_tb.append( s_v[i])
                else:                                    # Totally change table X values
                    for i in range(len(tb_x[par])):
                        par_tb.append(s_v[i])
                        par_tb.append(s_x[i])
                tmp_dict[par] = par_tb
                #print(tmp_dict[par])        
                theta_dict.update(tmp_dict)  
            else:
                raise Exception("There's something wrong with %s, please check it."%par)
        ##########################################################################
        cropdata.update(theta_dict)
        parameters = ParameterProvider(cropdata=cropdata,
                                       soildata=soil,
                                       sitedata=site)
        wofwof = Wofost71_PP(parameters,  weather, agromanagement)
        wofwof.run_till_terminate()    
        output = wofwof.get_output()
        summary_output = wofwof.get_summary_output()
        msg = "Reached maturity at {DOM} with max LAI of {LAIMAX} "\
    "and a yield of {TWSO} kg/ha."
        print(msg.format(**summary_output[0]))
        for var in varnames:
            tmp[var] = [t[var] for t in output]
        theta_dict["LAI"]=tmp["LAI"][-181:]
        theta_dict["day"]=tmp["day"][-181:]
        theta_dict["Yield"]=tmp["TWSO"][-1]
        outdata.append(theta_dict)
    np.save(out_en_file, outdata)

if __name__ == "__main__":

    ensemble_wofost(weather_type = "NASA", weather_path = "E:/temp/era5/", out_en_file = data_dir+"WOFOST_par_ensemble.npy")

