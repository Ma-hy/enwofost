use ECMWF ERA5 weather drive:

ensemble_wofost(lon = 115.55, lat=38., start = dt.date(2008,10,12),
                    end = None, en_size = 500, weather_type = "ERA", 
                    weather_path = %the folder you want to store the weahter data%, out_en_file = "WOFOST_par_ensemble.npy")

Ensemble of stuff. Wofost. Blah Blah 
You can now generate WOFOST ensemble in a give location quickly with 3 kinds of weather driver.
