#!/usr/bin/python
#=========================================================================
#!!!!!PLEASE CHANGE THIS NOTE IN CASE OF IMPROVING/CHANGING THIS CODE!!!!!
#VPRM AND ERA5 FOR T2M (0.1 DEG); NO NEED FOR REVERSE LATITUTE IN T2M
#=========================================================================
import numpy as np
from mpl_toolkits.basemap import Basemap
import math 
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
import netCDF4 as nc
from netCDF4 import Dataset,num2date, date2num
from datetime import datetime, timedelta
import datetime as dt
import sys 
STATIONS=['NOBS','HARV','HOW','LUCKY','TONZI','VAIRA','EAST','MEAD','OTHER']
main_dir='/home/georgy/WRF-CHEM/VPRM/'
Year='2019'
Month='03'
Date=Month+'-'+Year
DS_PAR = Dataset(main_dir+'Meteo/St_Petersburg/SWDOWN_d01_'+Month+'_'+Year+'_cut_mole_v2.nc')
DS_LSWI_1 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_14_cut_Int_v2.nc')
DS_LSWI_2 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_22_cut_Int_v2.nc')
DS_LSWI_3 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_30_cut_Int_v2.nc')
DS_LSWI_4 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_04_07_cut_Int_v2.nc')
DS_LSWI_5 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_04_15_cut_Int_v2.nc')
DS_LSWI_6 = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_04_23_cut_Int_v2.nc')
DS_VegType = Dataset(main_dir+'MODIS/MCD12Q1/MCD12Q1_Veg_Type_cut_Int_v2.nc')
i1=0
for d in range(1,32):
    if d<10:
        Day='0'+str(d)
        for h in range(0,24,3):
            if h<10:
                Hour='0'+str(h)
            else:
                Hour=str(h)
            j=0 
            NEE_MASS_lon=[]
            NEE_MASS=[]
            NEE_LAT=[]
            NEE_LON=[]   
            DS_T2M = Dataset(main_dir+'Meteo/St_Petersburg/T2M/Mar_2019/ERA5/ERA5_T2M_01_31-'+Date+'_025Deg_Int.nc')  
            if d>=1 and d<18:
                DS_LSWI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_14_cut_Int_v2.nc')
                DS_EVI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/EVI_vf/EVI_2019_03_14_cut_Int_v2.nc')                
            #lat_t2m=DS_T2M.variables['latitude'][:]
            lon=DS_PAR.variables['longitude'][:]
            lat=DS_PAR.variables['latitude'][:]
            T2M=DS_T2M.variables['t2m'][:]
            T2M=T2M-273.15
            T2M_TIME=DS_T2M.variables['time']
            T2M_ACTUAL_TIME_CONV=nc.num2date(T2M_TIME[:],T2M_TIME.units,T2M_TIME.calendar)
            PAR=DS_PAR.variables['SWDOWN'][:]
            LSWI=DS_LSWI.variables['Band1'][:]
            LSWI_1=DS_LSWI_1.variables['Band1'][:]
            LSWI_2=DS_LSWI_2.variables['Band1'][:]
            LSWI_3=DS_LSWI_3.variables['Band1'][:]
            LSWI_4=DS_LSWI_4.variables['Band1'][:]
            LSWI_5=DS_LSWI_5.variables['Band1'][:]
            LSWI_6=DS_LSWI_6.variables['Band1'][:]
            EVI=DS_EVI.variables['Band1'][:]
            VegType=DS_VegType.variables['Band1'][:]
            PAR_ACTUAL_DATE=DS_PAR.variables['XTIME']
            PAR_ACTUAL_DATE_CONV=nc.num2date(PAR_ACTUAL_DATE[:],PAR_ACTUAL_DATE.units,PAR_ACTUAL_DATE.calendar)
            i=0
            VT_TEST=[]
            for t in T2M_ACTUAL_TIME_CONV:
                if t.isoformat()[5:7]==Month and t.isoformat()[8:10]==Day and t.isoformat()[11:13]==Hour:
                    DATE=t
                    t1=i
                i=i+1    
            t=0
            i=0
            i2=0
            for t in PAR_ACTUAL_DATE_CONV:
                if t==DATE:
                    i1=i2
        #CONVERT ACCUMULATED PAR FROM um/m2 to um/m2s
                    if PAR_ACTUAL_DATE_CONV[i1].hour==3 or PAR_ACTUAL_DATE_CONV[i1].hour==15:
                        PAR_us=PAR[i1,:,:]#/(3*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==6 or PAR_ACTUAL_DATE_CONV[i1].hour==18:
                        PAR_us=PAR[i1,:,:]#/(6*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==9 or PAR_ACTUAL_DATE_CONV[i1].hour==21:
                        PAR_us=PAR[i1,:,:]#/(9*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==12 or PAR_ACTUAL_DATE_CONV[i1].hour==0:
                        PAR_us=PAR[i1,:,:]#/(12*3600)
                    t=0
                    for y in lat:
                        i=0
                        TEST=NEE_MASS_lon
                        NEE_MASS_lon=[]
                        for x in lon:
                            if math.isnan(PAR_us[j,i])==True:
                                PAR_US=0
                            else:
                                PAR_US=PAR_us[j,i]
                            VT=VegType[j,i]
                            VT_TEST.append(VT)
                            if VT==1 or VT==2:   #NOBS
                                VEGT=STATIONS[0] 
                                Tlow=1.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=262.0
                                alpha=0.244
                                beta=0.14
                                lambd=0.234
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                ##Pscale=(1+LSWI[j,i])/2
                                Pscale=1 #FOR EVERGREEN
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==3 or VT==4: #HARVARD
                                VEGT=STATIONS[1] 
                                Tlow=5.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=570.0
                                alpha=0.271
                                beta=0.25
                                lambd=0.127   
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==5:          #HOWLAND
                                VEGT=STATIONS[2]
                                Tlow=2.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=629.0
                                alpha=0.244
                                beta=-0.24
                                lambd=0.123 
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==6 or VT==7: #LUCKY HILLS
                                VEGT=STATIONS[3] 
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=321.0
                                alpha=0.028
                                beta=0.48
                                lambd=0.122
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==8 or VT==9: #TONZI RANCH
                                VEGT=STATIONS[4]
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=3241.0
                                alpha=0.012
                                beta=0.58
                                lambd=0.057
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==10:          #VAIRA
                                VEGT=STATIONS[5] 
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=18.0
                                PAR0=542.0
                                alpha=0.028
                                beta=0.72
                                lambd=0.213
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==11:          #EAST PEATLAND
                                VEGT=STATIONS[6]
                                Tlow=3
                                Tmin=0
                                Tmax=40
                                Topt=20
                                PAR0=558
                                alpha=0.081
                                beta=0.24
                                lambd=0.051 
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==12:            #MEAD
                                VEGT=STATIONS[7] 
                                Tlow=2.0
                                Tmin=5.0
                                Tmax=40.0
                                Topt=22.0
                                PAR0=11250.0
                                alpha=0.173
                                beta=0.82
                                lambd=0.075
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT>12 or np.ma.is_masked(VT)==True:             #OTHER
                                VEGT=STATIONS[8] 
                                NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            i=i+1
                        NEE_MASS_lon=np.array(NEE_MASS_lon)
                        #print(NEE_MASS_lon.shape)
                        NEE_MASS.append(NEE_MASS_lon)
                        j=j+1
                    #CREATE NETCDF WITH NEE
                    def newnc(lat,lon,Var,time): 
                        nrows=len(lat)
                        ncols=len(lon)
                        t=12
                        #t=1
                        print(nrows,ncols,t)
                        ppp =nc.Dataset('/home/georgy/WRF-CHEM/VPRM/NEE_res/Orig/Mar_2019/'+'VPRM_NEE_'+Day+Date+Hour+'.nc', 'w', format='NETCDF4')   
                        #DIMENSIONS
                        ppp.createDimension('time',t)   #1D dimension TIME
                        ppp.createDimension('lat', nrows)   #1D dimension LAT
                        ppp.createDimension('lon', ncols)   #1D dimension LON
                        #VARIABLES
                        longitude = ppp.createVariable('lon', 'f8', ('lon',))  #Longitude
                        longitude.standard_name = 'longitude'
                        longitude.long_name='longitude'
                        longitude.units = 'degrees_east'
                        longitude.axis = "X" 
                        
                        latitude = ppp.createVariable('lat' , 'f8', ('lat',))    #Latitude
                        latitude.standard_name = 'latitude'
                        latitude.long_name='latitude'
                        latitude.units = 'degrees_north'
                        latitude.axis = "Y"
                        
                        t=ppp.createVariable('time','f8',('time',)) #time
                        t.units = 'hours since 1900-01-01 00:00:00 UTC'
                        Var=np.array([Var,Var,Var,Var,Var,Var,Var,Var,Var,Var,Var,Var])
                        z = ppp.createVariable('NEE_CO2', 'f8', ('time','lat','lon')) #Values
                        z.standard_name = 'CO2_NEE'
                        z.long_name='NEE_flux_and_sink_VPRM'
                        z.units = 'micromol m-2 s-1'
                        z.coordinates="lon lat"
                        d=ppp.createVariable('date','i4',('time',)) #date
                        d.long_name="date"
                        d.units="YYYYMMDD"
                        ds=ppp.createVariable('datesec','i4',('time',)) #datesec
                        ds.long_name = "seconds_in_day"
                        ds[:]=[0,0,0,0,0,0,0,0,0,0,0,0]
                        d[:] = [20190101,20190201,20190301,20190401,20190501,20190601,20190701,20190801,20190901,20191001,20191101,20191201]
                        #WRITING VALUES TO THE NEW NETCDF VARIABLES
                        t[:]= time
                        latitude[:]= lat
                        longitude[:]= lon
                        #z[:,:] = np.array(Var)
                        z[:,:,:] = Var
                        ppp.close()   
                        return 0
                    print(newnc(lat,lon,NEE_MASS,PAR_ACTUAL_DATE[i1]))
                i2=i2+1
    else:
        Day=str(d)
        for h in range(0,24,3):
            if h<10:
                Hour='0'+str(h)
            else:
                Hour=str(h)     
            j=0 
            NEE_MASS_lon=[]
            NEE_MASS=[]
            NEE_LAT=[]
            NEE_LON=[]
            DS_T2M = Dataset(main_dir+'Meteo/St_Petersburg/T2M/Mar_2019/ERA5/ERA5_T2M_01_31-'+Date+'_025Deg_Int.nc')
            if d>=18 and d<26:
                DS_LSWI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_22_cut_Int_v2.nc')
                DS_EVI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/EVI_vf/EVI_2019_03_22_cut_Int_v2.nc')  
            elif d>=26:
                DS_LSWI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/LSWI_vf/LSWI_2019_03_30_cut_Int_v2.nc')
                DS_EVI = Dataset(main_dir+'MODIS/MOD09A1/St_Petersburg/EVI_vf/EVI_2019_03_30_cut_Int_v2.nc')              
            #lat_t2m=DS_T2M.variables['latitude'][:]
            lon=DS_PAR.variables['longitude'][:]
            lat=DS_PAR.variables['latitude'][:]
            T2M=DS_T2M.variables['t2m'][:]
            T2M=T2M-273.15
            T2M_TIME=DS_T2M.variables['time']
            T2M_ACTUAL_TIME_CONV=nc.num2date(T2M_TIME[:],T2M_TIME.units,T2M_TIME.calendar)
            PAR=DS_PAR.variables['SWDOWN'][:]
            LSWI=DS_LSWI.variables['Band1'][:]
            LSWI_1=DS_LSWI_1.variables['Band1'][:]
            LSWI_2=DS_LSWI_2.variables['Band1'][:]
            LSWI_3=DS_LSWI_3.variables['Band1'][:]
            LSWI_4=DS_LSWI_4.variables['Band1'][:]
            LSWI_5=DS_LSWI_5.variables['Band1'][:]
            LSWI_6=DS_LSWI_6.variables['Band1'][:]
            EVI=DS_EVI.variables['Band1'][:]
            VegType=DS_VegType.variables['Band1'][:]
            PAR_ACTUAL_DATE=DS_PAR.variables['XTIME']
            PAR_ACTUAL_DATE_CONV=nc.num2date(PAR_ACTUAL_DATE[:],PAR_ACTUAL_DATE.units,PAR_ACTUAL_DATE.calendar)
            i=0
            VT_TEST=[]
            for t in T2M_ACTUAL_TIME_CONV:
                if t.isoformat()[5:7]==Month and t.isoformat()[8:10]==Day and t.isoformat()[11:13]==Hour:
                    DATE=t
                    t1=i
                i=i+1    
            t=0
            i=0
            i2=0
            for t in PAR_ACTUAL_DATE_CONV:
                if t==DATE:
                    i1=i2
                    #CONVERT ACCUMULATED PAR FROM um/m2 to um/m2s
                    if PAR_ACTUAL_DATE_CONV[i1].hour==3 or PAR_ACTUAL_DATE_CONV[i1].hour==15:
                        PAR_us=PAR[i1,:,:]#/(3*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==6 or PAR_ACTUAL_DATE_CONV[i1].hour==18:
                        PAR_us=PAR[i1,:,:]#/(6*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==9 or PAR_ACTUAL_DATE_CONV[i1].hour==21:
                        PAR_us=PAR[i1,:,:]#/(9*3600)
                    elif PAR_ACTUAL_DATE_CONV[i1].hour==12 or PAR_ACTUAL_DATE_CONV[i1].hour==0:
                        PAR_us=PAR[i1,:,:]#/(12*3600)
                    t=0
                    for y in lat:
                        i=0
                        TEST=NEE_MASS_lon
                        NEE_MASS_lon=[]
                        for x in lon:
                            if math.isnan(PAR_us[j,i])==True:
                                PAR_US=0
                            else:
                                PAR_US=PAR_us[j,i]
                            VT=VegType[j,i]
                            VT_TEST.append(VT)
                            if VT==1 or VT==2:   #NOBS
                                VEGT=STATIONS[0] 
                                Tlow=1.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=262.0
                                alpha=0.244
                                beta=0.14
                                lambd=0.234
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                ##Pscale=(1+LSWI[j,i])/2
                                Pscale=1 #FOR EVERGREEN
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                                #=================================
                                if d==23 and j==13 and i==26:
                                    print('i1,t1,Tlow,Tmin,Tmax,Topt,T2M[t1,j,i],Tscale,Pscale,Wscale,alpha,beta,lambd,PAR0,LSWI_max,LSWI[j,i],EVI[j,i],PAR_US,GEE,R,NEE')
                                    print(i1,t1,Tlow,Tmin,Tmax,Topt,T2M[t1,13,26],Tscale,Pscale,Wscale,alpha,beta,lambd,PAR0,LSWI_max,LSWI[13,26],EVI[13,26],PAR_US,GEE,R,NEE)
                                #=================================
                            elif VT==3 or VT==4: #HARVARD
                                VEGT=STATIONS[1] 
                                Tlow=5.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=570.0
                                alpha=0.271
                                beta=0.25
                                lambd=0.127   
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==5:          #HOWLAND
                                VEGT=STATIONS[2]
                                Tlow=2.0
                                Tmin=0.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=629.0
                                alpha=0.244
                                beta=-0.24
                                lambd=0.123 
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==6 or VT==7: #LUCKY HILLS
                                VEGT=STATIONS[3] 
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=321.0
                                alpha=0.028
                                beta=0.48
                                lambd=0.122
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==8 or VT==9: #TONZI RANCH
                                VEGT=STATIONS[4]
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=20.0
                                PAR0=3241.0
                                alpha=0.012
                                beta=0.58
                                lambd=0.057
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==10:          #VAIRA
                                VEGT=STATIONS[5] 
                                Tlow=1.0
                                Tmin=2.0
                                Tmax=40.0
                                Topt=18.0
                                PAR0=542.0
                                alpha=0.028
                                beta=0.72
                                lambd=0.213
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==11:          #EAST PEATLAND
                                VEGT=STATIONS[6]
                                Tlow=3
                                Tmin=0
                                Tmax=40
                                Topt=20
                                PAR0=558
                                alpha=0.081
                                beta=0.24
                                lambd=0.051 
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT==12:            #MEAD
                                VEGT=STATIONS[7] 
                                Tlow=2.0
                                Tmin=5.0
                                Tmax=40.0
                                Topt=22.0
                                PAR0=11250.0
                                alpha=0.173
                                beta=0.82
                                lambd=0.075
                                #Tscale
                                if T2M[t1,j,i]>=Tmin:
                                    Tscale=((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax))/((T2M[t1,j,i]-Tmin)*(T2M[t1,j,i]-Tmax)-(T2M[t1,j,i]-Topt)**2)
                                else:
                                    Tscale=0
                                #Pscale
                                Pscale=(1+LSWI[j,i])/2
                                #Wscale
                                LSWI_max=max(LSWI_1[j,i],LSWI_2[j,i],LSWI_3[j,i],LSWI_4[j,i],LSWI_5[j,i],LSWI_6[j,i])
                                Wscale=(1+LSWI[j,i])/(1+LSWI_max)
                                #GEE
                                GEE=(lambd*Tscale*Pscale*Wscale*EVI[j,i]*PAR_US)/(1+(PAR_US/PAR0))
                                #R
                                R=alpha*T2M[t1,j,i]+beta
                                #NEE
                                NEE=-GEE+R
                                if math.isnan(LSWI[j,i])==True or math.isnan(LSWI_max)==True or math.isnan(EVI[j,i])==True:
                                    NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            elif VT>12 or np.ma.is_masked(VT)==True:             #OTHER
                                VEGT=STATIONS[8] 
                                NEE=0
                                #ADD TO ARRAY
                                NEE_MASS_lon.append(NEE)
                            i=i+1
                        NEE_MASS_lon=np.array(NEE_MASS_lon)
                        #print(NEE_MASS_lon.shape)
                        NEE_MASS.append(NEE_MASS_lon)
                        j=j+1      
                    #CREATE NETCDF WITH NEE
                    def newnc(lat,lon,Var,time): 
                        nrows=len(lat)
                        ncols=len(lon)
                        #t=1
                        t=12
                        print(nrows,ncols,t)
                        ppp =nc.Dataset('/home/georgy/WRF-CHEM/VPRM/NEE_res/Orig/Mar_2019/'+'VPRM_NEE_'+Day+Date+Hour+'.nc', 'w', format='NETCDF4')   
                        #DIMENSIONS
                        ppp.createDimension('time',t)   #1D dimension TIME
                        ppp.createDimension('lat', nrows)   #1D dimension LAT
                        ppp.createDimension('lon', ncols)   #1D dimension LON
                        #VARIABLES
                        longitude = ppp.createVariable('lon', 'f8', ('lon',))  #Longitude
                        longitude.standard_name = 'longitude'
                        longitude.long_name='longitude'
                        longitude.units = 'degrees_east'
                        longitude.axis = "X" 
                        
                        latitude = ppp.createVariable('lat' , 'f8', ('lat',))    #Latitude
                        latitude.standard_name = 'latitude'
                        latitude.long_name='latitude'
                        latitude.units = 'degrees_north'
                        latitude.axis = "Y"
                        
                        t=ppp.createVariable('time','f8',('time',)) #time
                        t.units = 'hours since 1900-01-01 00:00:00 UTC'
                        Var=np.array([Var,Var,Var,Var,Var,Var,Var,Var,Var,Var,Var,Var])
                        z = ppp.createVariable('NEE_CO2', 'f8', ('time','lat','lon')) #Values
                        z.standard_name = 'CO2_NEE'
                        z.long_name='NEE_flux_and_sink_VPRM'
                        z.units = 'micromol m-2 s-1'
                        z.coordinates="lon lat"
                        d=ppp.createVariable('date','i4',('time',)) #date
                        d.long_name="date"
                        d.units="YYYYMMDD"
                        ds=ppp.createVariable('datesec','i4',('time',)) #datesec
                        ds.long_name = "seconds_in_day"
                        ds[:]=[0,0,0,0,0,0,0,0,0,0,0,0]
                        d[:] = [20190101,20190201,20190301,20190401,20190501,20190601,20190701,20190801,20190901,20191001,20191101,20191201]
                        #WRITING VALUES TO THE NEW NETCDF VARIABLES
                        t[:]= time
                        latitude[:]= lat
                        longitude[:]= lon
                        #z[:,:] = np.array(Var)
                        z[:,:,:] = Var
                        ppp.close()   
                        return 0
                    print(newnc(lat,lon,NEE_MASS,PAR_ACTUAL_DATE[i1]))     
                i2=i2+1                        
