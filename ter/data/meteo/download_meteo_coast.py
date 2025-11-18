#globals().clear()


from pathlib import Path
import pandas as pd
from numpy import sqrt, mean
import datetime
import matplotlib.pyplot as plt
import requests
from datetime import datetime
import os
from numpy import sqrt, mean
import openmeteo_requests
import requests_cache
import pandas as pd
from retry_requests import retry


#Streamflow data from ACA
case_study = "ter"
homedir ="/home/dmercado/Documents/NEREIDA/"
os.chdir(homedir)

# historical precipitation data starts here
#set start date according to available data of fDOM, 
#because we will use daily data anyways
start_date = '1993-01-01'
start_date = pd.Timestamp(start_date)
end_date = '2025-06-30' #las day is not complete
end_date = pd.Timestamp(end_date)

#test = obs_data.resample('d').mean()

lat = 42.02
lon = 3.19

#Meteorological data
# Setup the Open-Meteo API client with cache and retry on error
cache_session = requests_cache.CachedSession('.cache', expire_after = -1)
retry_session = retry(cache_session, retries = 5, backoff_factor = 0.2)
openmeteo = openmeteo_requests.Client(session = retry_session)

# Make sure all required weather variables are listed here
# The order of variables in hourly or daily is important to assign them correctly below
url = "https://archive-api.open-meteo.com/v1/archive"
params = {
	"latitude": lat,
	"longitude": lon,
	"start_date": start_date.date().strftime('%Y-%m-%d'),
	"end_date": end_date.date().strftime('%Y-%m-%d'),
	"hourly": ["temperature_2m", "relative_humidity_2m",
	"precipitation", "surface_pressure", "cloud_cover",
	"et0_fao_evapotranspiration", "wind_speed_10m", "shortwave_radiation", 
	"soil_temperature_0_to_7cm", "soil_temperature_7_to_28cm", 
	"soil_temperature_28_to_100cm", "soil_temperature_100_to_255cm", 
	"soil_moisture_0_to_7cm", "soil_moisture_7_to_28cm", 
	"soil_moisture_28_to_100cm", "soil_moisture_100_to_255cm"],
	"models": "era5"
}
responses = openmeteo.weather_api(url, params=params)

# Process first location. Add a for-loop for multiple locations or weather models
response = responses[0]
print(f"Coordinates {response.Latitude()}°N {response.Longitude()}°E")
print(f"Elevation {response.Elevation()} m asl")
print(f"Timezone {response.Timezone()} {response.TimezoneAbbreviation()}")
print(f"Timezone difference to GMT+0 {response.UtcOffsetSeconds()} s")

# Process hourly data. The order of variables needs to be the same as requested.
hourly = response.Hourly()
hourly_temperature_2m = hourly.Variables(0).ValuesAsNumpy()
hourly_relative_humidity_2m = hourly.Variables(1).ValuesAsNumpy()
hourly_precipitation = hourly.Variables(2).ValuesAsNumpy()
hourly_surface_pressure = hourly.Variables(3).ValuesAsNumpy()
hourly_cloud_cover = hourly.Variables(4).ValuesAsNumpy()
hourly_et0_fao_evapotranspiration = hourly.Variables(5).ValuesAsNumpy()
hourly_wind_speed_10m = hourly.Variables(6).ValuesAsNumpy()
hourly_shortwave_radiation = hourly.Variables(7).ValuesAsNumpy()
hourly_soil_temperature_0_to_7cm = hourly.Variables(8).ValuesAsNumpy()
hourly_soil_temperature_7_to_28cm = hourly.Variables(9).ValuesAsNumpy()
hourly_soil_temperature_28_to_100cm = hourly.Variables(10).ValuesAsNumpy()
hourly_soil_temperature_100_to_255cm = hourly.Variables(11).ValuesAsNumpy()
hourly_soil_moisture_0_to_7cm = hourly.Variables(12).ValuesAsNumpy()
hourly_soil_moisture_7_to_28cm = hourly.Variables(13).ValuesAsNumpy()
hourly_soil_moisture_28_to_100cm = hourly.Variables(14).ValuesAsNumpy()
hourly_soil_moisture_100_to_255cm = hourly.Variables(15).ValuesAsNumpy()

hourly_data = {"date": pd.date_range(
	start = pd.to_datetime(hourly.Time(), unit = "s", utc = True),
	end = pd.to_datetime(hourly.TimeEnd(), unit = "s", utc = True),
	freq = pd.Timedelta(seconds = hourly.Interval()),
	inclusive = "left"
)}
hourly_data["temperature_2m"] = hourly_temperature_2m
hourly_data["relative_humidity_2m"] = hourly_relative_humidity_2m
hourly_data["precipitation"] = hourly_precipitation
hourly_data["surface_pressure"] = hourly_surface_pressure
hourly_data["cloud_cover"] = hourly_cloud_cover
hourly_data["et0_fao_evapotranspiration"] = hourly_et0_fao_evapotranspiration
hourly_data["wind_speed_10m"] = hourly_wind_speed_10m
hourly_data["shortwave_radiation"] = hourly_shortwave_radiation
hourly_data["soil_temperature_0_to_7cm"] = hourly_soil_temperature_0_to_7cm
hourly_data["soil_temperature_7_to_28cm"] = hourly_soil_temperature_7_to_28cm
hourly_data["soil_temperature_28_to_100cm"] = hourly_soil_temperature_28_to_100cm
hourly_data["soil_temperature_100_to_255cm"] = hourly_soil_temperature_100_to_255cm
hourly_data["soil_moisture_0_to_7cm"] = hourly_soil_moisture_0_to_7cm
hourly_data["soil_moisture_7_to_28cm"] = hourly_soil_moisture_7_to_28cm
hourly_data["soil_moisture_28_to_100cm"] = hourly_soil_moisture_28_to_100cm
hourly_data["soil_moisture_100_to_255cm"] = hourly_soil_moisture_100_to_255cm

hourly_dataframe = pd.DataFrame(data = hourly_data)
print(hourly_dataframe)

hourly_dataframe.to_csv(case_study+'/data/meteo/meteo-soil_data_coast.csv', index=False)

