# largely based on:
#  https://github.com/fvalka/nwp-sounding
#  https://github.com/fvalka/icon-skewt-plot

from metpy.units import units
import metpy.calc as mpcalc
import xarray as xr
import pytz
import bz2
import datetime
import json
import logging
import os
import os.path
import re
import numpy as np
import sys

import pint
ureg = pint.UnitRegistry()

# https://github.com/hgrecco/pint/issues/185
ureg = pint.UnitRegistry(preprocessors=[
    lambda s: s.replace('%%', ' permille '),
    lambda s: s.replace('%', ' percent '),
])
ureg.define(pint.unit.UnitDefinition(
    'permille', '%%', (), pint.converters.ScaleConverter(0.001),
))
ureg.define(pint.unit.UnitDefinition(
    'percent', '%', (), pint.converters.ScaleConverter(0.01),
))

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def level_path(prefix, model, scope, grid, run_hour, run_datetime, timestep, parameter, level, level_type):
    if level_type not in ("model-level", "time-invariant", "single-level"):
        raise AttributeError("Invalid level type")
    path = f"./{model}" \
        f"/grib" \
        f"/{run_hour:02d}" \
        f"/{parameter.lower()}" \
        f"/{model}_{scope}_{grid}_{level_type}_{run_datetime}_{timestep:03d}_{level}_{parameter.lower()}.grib2"
    return os.path.join(prefix, path)


class AllLevelDataResult:
    def __init__(self, data, model_time, valid_time):
        self.data = data.tolist()
        self.model_time = str(np.datetime_as_string(model_time))
        self.valid_time = str(np.datetime_as_string(valid_time))


def parameter_all_levels(prefix, model, scope, grid,
                         latitude, longitude, run_hour, run_datetime, timestep,
                         parameter, level_type="model-level",
                         base_level=60, top_level=1):

    levels = list(range(base_level, top_level - 1, -1))
    paths = [level_path(prefix, model, scope, grid, run_hour, run_datetime, timestep,
                        parameter, level, level_type) for level in levels]

    data_set = xr.open_mfdataset(paths,
                                 engine="cfgrib",
                                 concat_dim="generalVerticalLayer",
                                 combine='nested')

    interpolated = data_set.to_array()[0].interp(latitude=latitude, longitude=longitude)

    data = AllLevelDataResult(
        interpolated.values, interpolated.time.values, interpolated.valid_time.values)
    data_set.close()
    return data


def sounding(prefix, model, scope, grid, latitude, longitude, run_hour, run_datetime, timestep, parameter, level_type, base_level, top_level):

    sounding = parameter_all_levels(prefix, model, scope, grid, latitude, longitude,
                                    run_hour, run_datetime, int(timestep),
                                    parameter, level_type, base_level, top_level)
    response = sounding.__dict__
    return response


if __name__ == "__main__":

    data = {}
    params = ('t', 'p', 'qv', 'u', 'v')

    top = "/Users/mah/Ballon/src/downloader/2020121212"

    lat = 47.1026716552433
    lon = 15.2189647240443
    base_level = 65
    top_level = 1
    model = "icon-d2"
    grid = "regular-lat-lon"
    data['hhl'] = sounding(top, model, "germany", grid, lat, lon, 12, '2020121212', 0, 'hhl',
                           "time-invariant", base_level, top_level)

    for p in params:
        data[p] = sounding(top, model, "germany", grid, lat, lon, 12, '2020121212', 18, p,
                           "model-level", base_level, top_level)

    model_time = data['hhl']['model_time']
    valid_time = data['hhl']['valid_time']

    print(
        f"model={model}\nmodel_time={model_time}\nvalid_time={valid_time}\nlat={lat} lon={lon}\n")

    print(f"hhl\tHoehe\t\tTemp\t\tTau\t\t\trel Hum\t\t\tDruck\t\t\t\t\tWind\t\t\tWind aus")

    for i in range(top_level - 1, base_level - 1):
        t = data['t']['data'][i] * units.K
        p = (data['p']['data'][i]) * units.Pa
        qv = data['qv']['data'][i] * units("kg/kg")
        h = data['hhl']['data'][i] * units.meter
        u = data['u']['data'][i] * units('m/s')
        v = data['v']['data'][i] * units('m/s')

        # Dewpoint K
        dewpt = mpcalc.dewpoint_from_specific_humidity(qv, t, p)

        # relative humidity
        relhum = mpcalc.relative_humidity_from_dewpoint(t, dewpt)

        celsius = t.to('degC')
        ms = mpcalc.wind_speed(u, v)
        wdir = mpcalc.wind_direction(u, v, convention='from')

        dp = p.to('hPa')

        print(f"{i:-2d}\t{h:~5.1f}\t{celsius:~5.1f}\t{dewpt:~5.1f}\t{relhum:~%}\t{dp:6.1f}\t{ms:~3.1f}\t{wdir:3.0f}")
