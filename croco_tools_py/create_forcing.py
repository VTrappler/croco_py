#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import netCDF4 as netcdf
import time


def create_forcing(frcname, grdname, title,
                   smst, shft, swft, srft, sstt, ssst, smsc,
                   shfc, swfc, srfc, sstc, sssc):
    """
    Create an empty netcdf forcing file
    frcname: name of the forcing file
    grdname: name of the grid file
    title: title in the netcdf file
    """


    nc = netcdf.Dataset(grdname, 'r', format='NETCDF4')
    L = len(nc['xi_psi'])
    M = len(nc['eta_psi'])
    nc.close()
    Lp = L + 1
    Mp = M + 1

    nw = netcdf.Dataset(frcname, 'w', format='NETCDF4')

    xi_u = nw.createDimension('xi_u', L)
    eta_u = nw.createDimension('eta_u', Mp)
    xi_v = nw.createDimension('xi_v', Lp)
    eta_v = nw.createDimension('eta_v', M)
    xi_rho = nw.createDimension('xi_rho', Lp)
    eta_rho = nw.createDimension('eta_rho', Mp)
    xi_psi = nw.createDimension('xi_psi', L)
    eta_psi = nw.createDimension('eta_psi', M)
    sms_time = nw.createDimension('sms_time', len(smst))
    shf_time = nw.createDimension('shf_time', len(shft))
    swf_time = nw.createDimension('swf_time', len(swft))
    sst_time = nw.createDimension('sst_time', len(sstt))
    srf_time = nw.createDimension('srf_time', len(srft))
    sss_time = nw.createDimension('sss_time', len(ssst))
    wwv_time = nw.createDimension('wwv_time', len(wwvt))

    sms_time = nw.createVariable('sms_time', 'f8')
    sms_time.long_name = 'surface momentum stress time'
    sms_time.units = 'days'
    sms_time.cycle_length = smsc

    shf_time = nw.createVariable('shf_time', 'f8')
    shf_time.long_name = 'surface heat flux time'
    shf_time.units = 'days'
    shf_time.cycle_length = shfc 

    swf_time = nw.createVariable('swf_time', 'f8')
    swf_time.long_name = 'surface freshwater flux time'
    swf_time.units = 'days'
    swf_time.cycle_length = swfc

    sst_time = nw.createVariable('sst_time', 'f8')
    sst_time.long_name = 'sea surface temperature time'
    sst_time.units = 'days'
    sst_time.cycle_length = sstc

    sss_time = nw.createVariable('sss_time', 'f8')
    sss_time.long_name = 'sea surface salinity time'
    sss_time.units = 'days'
    sss_time.cycle_length = sssc

    srf_time = nw.createVariable('srf_time', 'f8')
    srf_time.long_name = 'solar shortwave radiation time'
    srf_time.units = 'days'
    srf_time.cycle_length = srfc

    wwv_time = nw.createVariable('wwv_time', 'f8')
    wwv_time.long_name = 'surface wave fields time'
    wwv_time.units = 'days'
    wwv_time.cycle_length = smsc

    sustr = nw.createVariable('sms_time', 'eta_u', 'xi_u', 'f8')
    sustr.long_name = 'surface u-momentum stress'
    sustr.units = 'Newton meter-2'

    svstr = nw.createVariable('sms_time', 'eta_v', 'xi_v', 'f8')
    svstr.long_name = 'surface v-momentum stress'
    svstr.units = 'Newton meter-2'

    shflux = nw.createVariable('shf_time', 'eta_rho', 'xi_rho', 'f8')
    shflux.long_name = 'surface net heat flux'
    shflux.units = 'Watts meter-2'

    swflux = nw.createVariable('swf_time', 'eta_rho', 'xi_rho', 'f8')
    swflux.long_name = 'surface freshwater flux (E-P)'
    swflux.units = 'centimeter day-1'
    swflux.positive = 'net evaporation'
    swflux.negative = 'net precipitation'

    SST = nw.createVariable('sst_time', 'eta_rho', 'xi_rho', 'f8')
    SST.long_name = 'sea surface temperature'
    SST.units = 'Celsius'

    SSS = nw.createVariable('sss_time', 'eta_rho', 'xi_rho', 'f8')
    SSS.long_name = 'sea surface salinity'
    SSS.units = 'PSU'

    dQdSST = nw.createVariable('sst_time', 'eta_rho', 'xi_rho', 'f8')
    dQdSST.long_name = 'surface net heat flux sensitivity to SST'
    dQdSST.units = 'Watts meter-2 Celsius-1'

    swrad = nw.createVariable('srf_time', 'eta_rho', 'xi_rho', 'f8')
    swrad.long_name = 'solar shortwave radiation'
    swrad.units = 'Watts meter-2'
    swrad.positive = 'downward flux, heating'
    swrad.negative = 'upward flux, cooling'

    Awave = nw.createVariable('wwv_time', 'eta_rho', 'xi_rho', 'f8')
    Awave.long_name = 'wind induced wave amplitude'
    Awave.units = 'm'

    Dwave = nw.createVariable('wwv_time', 'eta_rho', 'xi_rho', 'f8')
    Dwave.long_name = 'wind induced wave direction'
    Dwave.units = 'degree'

    Pwave = nw.createVariable('wwv_time', 'eta_rho', 'xi_rho', 'f8')
    Pwave.long_name = 'wind induced wave period'
    Pwave.units = 'second'
