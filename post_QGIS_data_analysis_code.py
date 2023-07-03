# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:22:09 2023

@author: mjrobinson
"""
#%%


# below is an example code of 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math




from matplotlib import rcParams

# from scipy.optimize import curve_fit
rcParams['font.family']='sans-serif'
rcParams['font.sans-serif']=['Arial']

# import functools as fnc
plt.rcParams.update({'font.size': 10})
palette = sns.color_palette("inferno", 9)
ccolor=[palette[3],palette[5],palette[6]]




#%%
# import data

filelist=["gabilan/","OCR/","WestVirgina/","san_gab/","Kauai/","utah3/"]
namelist=["Gabilan","OCR","West Virgina","San Gab","Kauai","Utah"]

numberlist=[1,2,3,4,5,6]
gcavg=np.zeros(6)
gcstd=np.zeros(6)
minareaperlandscae=np.zeros( (np.size(filelist)))
maxareaperlandscae=np.zeros( (np.size(filelist)))
mingc=np.zeros( (np.size(filelist)))
maxgc=np.zeros( (np.size(filelist)))


ocr_catagory=np.loadtxt("OCR/catagory.txt")

uplift=(.147,.1,2,.3,.03,7,.055,.5,)    
stdk=np.zeros(np.size(filelist))
gcmean=np.zeros(np.size(namelist))
gcstd=np.zeros(np.size(namelist))
gcfull=np.array(())

for h in range(np.size(filelist)):

    
  
    internalcount=1
    df=pd.read_csv('{}basins.csv'.format(filelist[h]),header=0,float_precision='round_trip')
    
    df=df.rename(columns={'fid': "basin_id"})
    basin_data = df
    del basin_data['DN']
    del basin_data['vertex_index']
    del basin_data['vertex_part']
    del basin_data['vertex_part_index']
    del basin_data['distance']
    del basin_data['angle']
    
    
    basin_data.columns =['basin_id','x','y','z']
    basin_data= basin_data.dropna()
    
    df1=pd.read_csv('{}channel1.csv'.format(filelist[h]),header=0,float_precision='round_trip')

    df1=df1.rename(columns={'fid': "basin_id"})
    df1=df1.rename(columns={'elevation': "z"})
    channel_data =df1
    del channel_data['feature_x']
    del channel_data['feature_y']
    del channel_data['nearest_x']
    del channel_data['nearest_y']
    del channel_data['distance']
    del channel_data['basin_key']
    del channel_data['latitude']
    del channel_data['longitude']
    del channel_data['StreamOrder']
    del channel_data['n']
    del channel_data['DN']
    
    
    channel_data.columns =['chi','drainage_area','z','flowdistance','source','basin_id','x','y']
    
    #
    
    
    channel_data_round = channel_data.astype('int')
    basin_data_round = basin_data.astype('int')
    
    
    #randomly generate basins to be analyzed this works but go with 25 biggest
    
    
    
    rng=channel_data_round.loc[channel_data_round.groupby('basin_id')['drainage_area'].nlargest(1).reset_index(0).index]
    # rng = rng.nlargest(20,['drainage_area'])
    basin_list=np.array(rng['basin_id'])
    print("this is {} basin list".format(namelist[h]))

       
    #



    #save longest flow distance for later and create some arrays to save stuff to 
    
    saveflowdist=np.zeros( (1,np.size(basin_list)+1))
    saveperimeter=np.zeros( np.size(basin_list))
    savearea=np.zeros( np.size(basin_list))
    axisratio=np.zeros((np.size(basin_list))*2)
    axisratiosingle=np.zeros((np.size(basin_list)))
    centery =  np.zeros((np.size(basin_list))) 
    gcthisbasin=np.zeros( np.size(basin_list))




    
    # loop through basin list for each landscape pulling out some iniatial area and  GC information.
  
    for f in range(0,np.size(basin_list),1):
        
    
       
        basinpick= basin_list[f]
        

        

                #stuff for perimeter
        for key, grp in channel_data.groupby(['basin_id']):
            if key == basinpick:
                savearea[f] =grp["drainage_area"].max()
        for key, grp2 in basin_data.groupby(['basin_id']):
            if key == basinpick:
                fakexperim = np.array(grp2['x'])
                fakeyperim = np.array(grp2['y'])
                around=0
                fakexridge = np.array(grp2['x'])
                fakeyridge = np.array(grp2['y'])
                fakezridge = np.array(grp2['z'])
                for i in range(np.size(fakexperim)):
                    if i == np.size(fakexperim)-1:
                        continue
                    around = around+np.sqrt(((fakexperim[i+1])-(fakexperim[i]))**2+((fakeyperim[i+1])-(fakeyperim[i]))**2)
    
  
                
        
        saveperimeter[f]=around
        gc=saveperimeter[f]/(2*np.sqrt(np.pi*savearea[f]))
        gcfull=np.append(gcfull,gc)
        gcthisbasin[f]=gc
        


  

            
    # checks
    minareaperlandscae[h]=np.min(savearea)
    maxareaperlandscae[h]=np.max(savearea)
    check=gcthisbasin[gcthisbasin>1]
    mingc[h]=np.min(check)
    maxgc[h]=np.max(check)
        
    #
    if h==0:
        gc_gabilan=saveperimeter/(2*np.sqrt(np.pi*savearea))
       
        basin_gabilan=basin_list
        gcavg[h]=np.mean(gc_gabilan)
        gcstd[h]=np.std(gc_gabilan)
    if h==1:

        gc_OCR=saveperimeter/(2*np.sqrt(np.pi*savearea))
        
        basin_OCR=basin_list
        gcavg[h]=np.mean(gc_OCR)
        gcstd[h]=np.std(gc_OCR)
        
    if h==2:
        gc_WestVirgina=saveperimeter/(2*np.sqrt(np.pi*savearea))
     
        basin_WestVirgina=basin_list
        gcavg[h]=np.mean(gc_WestVirgina)
        gcstd[h]=np.std(gc_WestVirgina)
        
        
    if h ==3:
        
        
        gc_san_gab=saveperimeter/(2*np.sqrt(np.pi*savearea))
        
        basin_san_gab=basin_list
        gcavg[h]=np.mean(gc_san_gab)
        gcstd[h]=np.std(gc_san_gab)
    if h ==4:
        
        
        gc_kauai=saveperimeter/(2*np.sqrt(np.pi*savearea))
      
        basin_kauai=basin_list
        gcavg[h]=np.mean(gc_kauai)
        gcstd[h]=np.std(gc_kauai)
    if h ==5:
        
        
        gc_utah=saveperimeter/(2*np.sqrt(np.pi*savearea))
        ar_utah=axisratiosingle
        basin_utah=basin_list
        gcavg[h]=np.mean(gc_utah)
        gcstd[h]=np.std(gc_utah)
     

      # create list for saving GC info 
            
gcnamelist=["gc_gabilan","gc_OCR","gc_WestVirgina","gc_san_gab","gc_kauai","gc_utah"]

        # save some inforation that can be used as checks later
   
   
np.savetxt("minarea.csv",minareaperlandscae)
np.savetxt("maxarea.csv",maxareaperlandscae)
np.savetxt("mingc.csv",mingc)
np.savetxt("maxgc.csv",maxgc)



#  I hardcoded a bunch of this to save time initially and never came back to polish it off
# it works but can be simplifed greatly. 

savealltribsx=np.array(())

savealltribsy=np.array(())
savealltriborder=np.array(())
savealltribgc=np.array(())
savealltribar=np.array(())
savealltribbasin=np.array(())
saveallname=np.array(())

basintotalarea=np.array(())
basindn1area=np.array(())
basindn2area=np.array(())
basindn3area=np.array(())

s1savealltribsx=np.array(())

s1savealltribsy=np.array(())
s1savealltriborder=np.array(())
s1savealltribgc=np.array(())
s1savealltribar=np.array(())
s1savealltribbasin=np.array(())
s1saveallname=np.array(())
s1savealltribsx=np.array(())
s1savetribarea=np.array(())

s2savealltribsx=np.array(())

s2savealltribsy=np.array(())
s2savealltriborder=np.array(())
s2savealltribgc=np.array(())
s2savealltribar=np.array(())
s2savealltribbasin=np.array(())
s2saveallname=np.array(())
s2savealltribsx=np.array(())
s2savetribarea=np.array(())


s3savealltribsx=np.array(())

s3savealltribsy=np.array(())
s3savealltriborder=np.array(())
s3savealltribgc=np.array(())
s3savealltribar=np.array(())
s3savealltribbasin=np.array(())
s3saveallname=np.array(())
s3savealltribsx=np.array(())
savetribarea=np.array(())
s3savetribarea=np.array(())

space1=np.array(())
space2=np.array(())
space3=np.array(())

slope1=np.array(())
slope2=np.array(())
slope3=np.array(())
trib0=np.array(())
trib1=np.array(())
trib2=np.array(())
trib3=np.array(())
trib0x=np.array(())
trib1x=np.array(())
trib2x=np.array(())
trib3x=np.array(())
trib0y=np.array(())
trib1y=np.array(())
trib2y=np.array(())
trib3y=np.array(())

newy=np.zeros(20)
newx=np.zeros(20)
newstd=np.zeros(20)


# now loop through basins for each landscape and pull out the meat of the inforation we want. trib spacing length and minstem info

for h in range(np.size(namelist)):
    area=np.array(())
    length=np.array(())

    if h==0:
        
           
        # df=pd.read_csv('{}try2.csv'.format(filelist[h]),header=0,float_precision='round_trip')
      
        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
     
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        # del df['StreamOrder']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})

        gcnew=gc_gabilan
 
            
        basin_list=basin_gabilan
        basin_list=basin_list[np.where(gcnew > 1)]
   
          
        gcnew=gcnew[np.where(gcnew> 1)]
        #work on only taking certain basins 

        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()

        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5
    if h==1:
        # df=pd.read_csv('{}try2.csv'.format(filelist[h]),header=0,float_precision='round_trip')
        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
       
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        # del df['StreamOrder']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})
        
        # df=df.rename(columns={'Stream Order': "Stream_Order"})
        gcnew=gc_OCR
     
        basin_list=basin_OCR
        basin_list=basin_list[np.where(gcnew > 1)]
       
        gcnew=gcnew[np.where(gcnew> 1)]
        
        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()

        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5
    if h ==2:

        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
       
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        # del df['StreamOrder']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})
        gcnew=gc_WestVirgina
      
        basin_list= basin_WestVirgina
        basin_list=basin_list[np.where(gcnew > 1)]
       
        gcnew=gcnew[np.where(gcnew> 1)]
        
        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()

        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5
       
    if h==3:

        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
     
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        # del df['StreamOrder']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})
          
        gcnew=gc_san_gab
        basin_list=basin_san_gab
       
        basin_list=basin_list[np.where(gcnew > 1)]
      
        gcnew=gcnew[np.where(gcnew> 1)]
        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()
  
        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5
     
    if h ==4:
        #%%
        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
     
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        del df['DN']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})
          
        gcnew=gc_kauai
        basin_list=basin_kauai
       
        basin_list=basin_list[np.where(gcnew > 1)]
      
        gcnew=gcnew[np.where(gcnew> 1)]
        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()
  
        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5
        #
    if h ==5:

        df=pd.read_csv("{}channel1.csv".format(filelist[h]),header=0,float_precision='round_trip')
        df=df.rename(columns={'fid': "FID_basin"})
        df=df.rename(columns={'elevation': "z"})
     
        del df['feature_x']
        del df['feature_y']
        del df['nearest_x']
        del df['nearest_y']
        del df['distance']
        del df['basin_key']
        del df['latitude']
        del df['longitude']
        del df['DN']
        del df['n']
        
        df=df.rename(columns={'StreamOrder': "Stream_Order"})


        gcnew=gc_utah
      
        basin_list=basin_utah
        basin_list=basin_list[np.where(gcnew > 1)]
   
        gcnew=gcnew[np.where(gcnew> 1)]
        #work on only taking certain basins 

        df=df[df["FID_basin"].isin(basin_list)]
        bybasin=df.groupby('FID_basin')['source_key'].aggregate(['min','max'])
        letssee=df.groupby("FID_basin")
        groups=df.groupby('FID_basin').groups
        basinid=groups.keys()
        basinid=list(basinid)
        basinid=np.array(basinid)
        how_many_basins=len(bybasin)
        tofind = df["FID_basin"].unique()

        
        minsource=np.array(bybasin['min'])
        maxsource=np.array(bybasin['max'])
        df["DSO"]=5

    
    
    
    



        
   
    np.savetxt("{}sourcelist.csv".format(namelist[h]),minsource)
    filtermainstem=[]
    ##the flow distance start at outlet=0 and goes to headwaters =max
    check1=0
    check2=0
    check3=0
  
    for i in range(np.size(basin_list)):
        newdn1areathing=np.array(())
        newdn2areathing=np.array(())
        newdn3areathing=np.array(())
        basin=basin_list[i]
        
        gc=gcnew[i]
 

                    #
        basinspecific1=np.array(())
        basinspecific2=np.array(())
        basinspecific3=np.array(())
        basinspecific1y=np.array(())
        basinspecific2y=np.array(())
        basinspecific3y=np.array(())
        sourcelist=np.arange((minsource[i]),maxsource[i]+1,1)
        tribgc =np.zeros(np.size(sourcelist))
        
        tribbasin=np.zeros(np.size(sourcelist))
        triblength = np.zeros(np.size(sourcelist))
        tribconfluencedistance=np.zeros(np.size(sourcelist))
        triborder=np.zeros(np.size(sourcelist))

        source1=df.loc[df['source_key'] == minsource[i]]
        mainstemx=np.array(source1["x"])
        mainstemy=np.array(source1["y"])
        

            #



        minus=source1.flow_distance.min()
        maxsorder=source1.Stream_Order.max()

        normalize=(source1.flow_distance.max())-minus
        if normalize == 0: 
            continue
        basintotalarea=np.append(basintotalarea,source1.drainage_area.max())
            
       
        
        for f in range(np.size(sourcelist)):
            
            trib=df.loc[df['source_key'] == sourcelist[f]]
     
            
            
            tribx=np.array(trib["x"][trib.flow_distance==np.min(trib.flow_distance)])
            triby=np.array(trib["y"][trib.flow_distance==np.min(trib.flow_distance)])

            #
            distance=[]
            
            for l in range(np.size(mainstemx)):
                fff=math.dist([tribx,triby],[mainstemx[l],mainstemy[l]])
                distance=np.append(distance,fff)
                #
            if np.min(distance)>100:
                triborder[f]=np.nan
                bottom=np.nan
                tribconfluencedistance[f]=np.nan
                triblength[f]=np.nan
                tribgc[f]=np.nan
          
                tribbasin[f]=np.nan
                
            

 
            
            else:
                filtermainstem=np.append(filtermainstem,sourcelist[f])
                triborder[f]=(maxsorder-trib.Stream_Order.max())
                bottom=trib.flow_distance.min()
                tribconfluencedistance[f]=((bottom- minus)/normalize)
                triblength[f]=(((trib.flow_distance.max())-bottom)/normalize)
                tribgc[f]=gcnew[i]
          
                tribbasin[f]=basin_list[i]
                
                df.DSO[df['source_key'] == sourcelist[f]]=triborder[f]
                
            if triborder[f]==1:
                newdn1areathing=np.append(newdn1areathing,trib.drainage_area.max())
            if triborder[f]==2:
                newdn2areathing=np.append(newdn2areathing,trib.drainage_area.max())
            if triborder[f]==3:
                newdn3areathing=np.append(newdn3areathing,trib.drainage_area.max())
        addareatolength=np.full(np.size(triblength),source1.drainage_area.max())
        nametoadd=np.full(np.size(triblength),numberlist[h])
        Idtoadd=np.full(np.size(triblength),numberlist[h])
    
        savealltribsx=np.append(savealltribsx,tribconfluencedistance)
        savealltribsy=np.append(savealltribsy,triblength)
        savealltriborder=np.append(savealltriborder,triborder)
        savetribarea=np.append(savetribarea,addareatolength)
        
        
        savealltribgc=np.append(savealltribgc,tribgc)
  
        savealltribbasin=np.append(savealltribbasin,tribbasin)
        saveallname=np.append(saveallname,nametoadd)
        #
        
        #SPACING STUFF first get stream order info
        
        x1=tribconfluencedistance[triborder==1]
        x1=1-x1
        
        savesorting1=x1
        x1=np.sort(x1)
        
        
        x2=tribconfluencedistance[triborder==2]
        x2=1-x2
        
        savesorting2=x2
        x2=np.sort(x2)
        
        
        x3=tribconfluencedistance[triborder==3]
        x3=1-x3
        savesorting3=x3
        x3=np.sort(x3)
        
        
        if np.size(x1)>2:
           tribspace1=x1[1:]-x1[:-1]
           space1=np.append(space1,tribspace1)
        if np.size(x1)==2:
         
            tribspace1=np.nan
            space1=np.append(space1,tribspace1)
        if np.size(x2)>2:
           tribspace2=x2[1:]-x2[:-1]
           space2=np.append(space2,tribspace2)
        if np.size(x2)==2:
         
            tribspace2=np.nan
            space2=np.append(space2,tribspace2)
            
            
        if np.size(x3)>2:
           tribspace3=x3[1:]-x3[:-1]
           space3=np.append(space3,tribspace3)
        if np.size(x3)==2:
          
            tribspace3=np.nan
            space3=np.append(space3,tribspace3)
   #
        
         # super bad coding but here it goes
        trib1area=addareatolength[triborder==1]
        trib1area=trib1area[np.argsort(savesorting1)]
         
        tribcon1=tribconfluencedistance[triborder==1]
        tribcon1=tribcon1[np.argsort(savesorting1)]
        
        tribL1=triblength[triborder==1]
        tribL1=tribL1[np.argsort(savesorting1)]
        
        tribeltaN1=triborder[triborder==1]
        tribeltaN1=tribeltaN1[np.argsort(savesorting1)]
        
        tribshape1=tribgc[triborder==1]
        tribshape1=tribshape1[np.argsort(savesorting1)]
        
        tribbasinidentify1=tribbasin[triborder==1]
        tribbasinidentify1=tribbasinidentify1[np.argsort(savesorting1)]
        
        triblandscapename1=nametoadd[triborder==1]
        triblandscapename1=triblandscapename1[np.argsort(savesorting1)]
        
        s1savealltribsx=np.append(s1savealltribsx,tribcon1[1:])
        s1savealltribsy=np.append(s1savealltribsy,tribL1[1:])
        s1savealltriborder=np.append(s1savealltriborder,tribeltaN1[1:])
        
        
        s1savealltribgc=np.append(s1savealltribgc,tribshape1[1:])
          
        s1savealltribbasin=np.append(s1savealltribbasin,tribbasinidentify1[1:])
        s1saveallname=np.append(s1saveallname,triblandscapename1[1:])
        s1savetribarea=np.append(s1savetribarea,trib1area[1:])
#
        trib2area=addareatolength[triborder==2]
        trib2area=trib2area[np.argsort(savesorting2)]

      
        tribcon2=tribconfluencedistance[triborder==2]
        tribcon2=tribcon2[np.argsort(savesorting2)]
        
        tribL2=triblength[triborder==2]
        tribL2=tribL2[np.argsort(savesorting2)]
        
        tribeltaN2=triborder[triborder==2]
        tribeltaN2=tribeltaN2[np.argsort(savesorting2)]
        
        tribshape2=tribgc[triborder==2]
        tribshape2=tribshape2[np.argsort(savesorting2)]
        
        tribbasinidentify2=tribbasin[triborder==2]
        tribbasinidentify2=tribbasinidentify2[np.argsort(savesorting2)]
        
        triblandscapename2=nametoadd[triborder==2]
        triblandscapename2=triblandscapename2[np.argsort(savesorting2)]
        
        s2savealltribsx=np.append(s2savealltribsx,tribcon2[1:])
        s2savealltribsy=np.append(s2savealltribsy,tribL2[1:])
        s2savealltriborder=np.append(s2savealltriborder,tribeltaN2[1:])
        
        
        s2savealltribgc=np.append(s2savealltribgc,tribshape2[1:])
          
        s2savealltribbasin=np.append(s2savealltribbasin,tribbasinidentify2[1:])
        s2saveallname=np.append(s2saveallname,triblandscapename2[1:])
        s2savetribarea=np.append(s2savetribarea,trib2area[1:])
#


        trib3area=addareatolength[triborder==3]
        trib3area=trib3area[np.argsort(savesorting3)]
        tribcon3=tribconfluencedistance[triborder==3]
        tribcon3=tribcon3[np.argsort(savesorting3)]
        
        tribL3=triblength[triborder==3]
        tribL3=tribL3[np.argsort(savesorting3)]
        
        tribeltaN3=triborder[triborder==3]
        tribeltaN3=tribeltaN3[np.argsort(savesorting3)]
        
        tribshape3=tribgc[triborder==3]
        tribshape3=tribshape3[np.argsort(savesorting3)]
        
        tribbasinidentify3=tribbasin[triborder==3]
        tribbasinidentify3=tribbasinidentify3[np.argsort(savesorting3)]
        
        triblandscapename3=nametoadd[triborder==3]
        triblandscapename3=triblandscapename3[np.argsort(savesorting3)]
        
        s3savealltribsx=np.append(s3savealltribsx,tribcon3[1:])
        s3savealltribsy=np.append(s3savealltribsy,tribL3[1:])
        s3savealltriborder=np.append(s3savealltriborder,tribeltaN3[1:])
        
        
        s3savealltribgc=np.append(s3savealltribgc,tribshape3[1:])
          
        s3savealltribbasin=np.append(s3savealltribbasin,tribbasinidentify3[1:])
        s3saveallname=np.append(s3saveallname,triblandscapename3[1:])
        s3savetribarea=np.append(s3savetribarea,trib3area[1:])        
           
  
        basindn1area=np.append(basindn1area,np.sum(newdn1areathing))     
        basindn2area=np.append(basindn2area,np.sum(newdn2areathing)) 
        basindn3area=np.append(basindn3area,np.sum(newdn3areathing)) 
        
        #
    df.to_csv("{}SAVETHISTHING.csv".format(namelist[h]))
    df2=df[df["source_key"].isin(filtermainstem)]
    df2.to_csv("{}mainstemfiltercheck.csv".format(namelist[h]))

#%%
#how much area is on average is made up of dn1 ect ....
tellmethis=np.mean(basindn1area/basintotalarea)
tellmethis2=np.mean(basindn2area/basintotalarea)
tellmethis3=np.mean(basindn3area/basintotalarea)


#%%
#more hard coding for spacing 

columnss=["basin_id","trib_length","along_distance","DSO","gc","spacing","landscape"]
s1=pd.DataFrame()
s1["along_distance"]=s1savealltribsx
s1["trib_length"]=s1savealltribsy
s1["DSO"]=s1savealltriborder
s1["landscape"]=s1saveallname
s1["gc"]=s1savealltribgc
s1["basin_id"]=s1savealltribbasin
s1["spacing"]=space1
s1["area"]=s1savetribarea
#%%
s2=pd.DataFrame(columns = columnss)
s2["along_distance"]=s2savealltribsx
s2["trib_length"]=s2savealltribsy
s2["DSO"]=s2savealltriborder
s2["landscape"]=s2saveallname
s2["gc"]=s2savealltribgc
s2["basin_id"]=s2savealltribbasin
s2["spacing"]=space2
s2["area"]=s2savetribarea
#%%

columnss=["basin_id","trib_length","along_distance","DSO","gc"]
tryy=pd.DataFrame(columns = columnss)
tryy["along_distance"]=savealltribsx
tryy["trib_length"]=savealltribsy
tryy["DSO"]=savealltriborder
tryy["landscape"]=saveallname
tryy["gc"]=savealltribgc
tryy["basin_id"]=savealltribbasin
tryy["area"]=savetribarea
tryy=tryy[tryy["DSO"]<=3]


s3=pd.DataFrame(columns = columnss)
s3["along_distance"]=s3savealltribsx
s3["trib_length"]=s3savealltribsy
s3["DSO"]=s3savealltriborder
s3["landscape"]=s3saveallname
s3["gc"]=s3savealltribgc
s3["basin_id"]=s3savealltribbasin
s3["spacing"]=space3
s3["area"]=s3savetribarea

#%%
print(len(tryy[(tryy["DSO"]==1)]))
tryy0=tryy[(tryy["DSO"]==0)]
tryy1=tryy[(tryy["DSO"]==1)]
tryy2=tryy[(tryy["DSO"]==2) ]
tryy3=tryy[(tryy["DSO"]==3) ]

#%%
 
tryy.to_csv("tryy.csv")
s1.to_csv("s1.csv")
s2.to_csv("s2.csv")
s3.to_csv("s3.csv")
#%%
tryy=pd.read_csv("tryy.csv")
s1=pd.read_csv("s1.csv")
s2=pd.read_csv("s2.csv")
s3=pd.read_csv("s3.csv")

print(len(tryy[(tryy["DSO"]==1)]))
tryy0=tryy[(tryy["DSO"]==0)]
tryy1=tryy[(tryy["DSO"]==1)]
tryy2=tryy[(tryy["DSO"]==2) ]
tryy3=tryy[(tryy["DSO"]==3) ]






