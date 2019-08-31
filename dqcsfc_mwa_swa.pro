PRO dqcsfc_mwa_swa
  COMPILE_OPT IDL2
  e=ENVI(/HEADLESS);开始一个无窗口的应用程序
  ;第一部分:辐射亮度定标.打开头文件
  MTLFile_address=ENVI_PICKFILE(TITLE='请选择Landsat 8 MTL文件',FILTER='*_MTL.txt')
  Dirname=FILE_DIRNAME(MTLFile_address)
  RawRaster=e.OpenRaster(MTLFile_address);通过地址打开这个数据集,一共五部分数据集
  Mul_RawRaster=RawRaster[0];获得第一部分数据集,也就是多光谱部分,共7个波段
  MulRaw_fid=ENVIRasterToFID(Mul_RawRaster)
  Tir_RawRaster=RawRaster[3];获得第四部分数据集,也就是热红外部分,共2个波段
  TirRaw_fid=ENVIRasterToFID(Tir_RawRaster)
  MTLfile_name=FILE_BASENAME(MTLFile_address);通过地址得到数据文件名
  pPos=STRPOS(MTLfile_name,'.',/reverse_search);pPos为'.'的位置
  ;imageprefix为原先去除后缀后的文件名,去掉了MLT.dat和其他后缀
  imageprefix=strmid(MTLfile_name,0,pPos-18)
  ;获取现在的儒略历时间,换算为公历,将时间存入内存
  caldat,systime(/JULIAN),month,day,year,hour,min,sec
  year=strtrim(year,2);将年份去掉前后空格
  if((month ge 1)and(month le 9))then begin;如果月份在1月到9月之间
    month=strjoin(['0',strtrim(month,2)]);把月份前后的空格去掉,在前面加0
  endif else begin
    month=strtrim(month,2);否则无需加0,只把月份前后的空格去掉即可
  endelse
  if((day ge 1)and(day le 9))then begin;如果日期在1日到9日之间
    day=strjoin(['0',strtrim(day,2)])
  endif else begin
    day=strtrim(day,2)
  endelse
  if((hour ge 0)and(hour le 9))then begin;如果小时在0时到9时之间
    hour=strjoin(['0',strtrim(hour,2)])
  endif else begin
    hour=strtrim(hour,2)
  endelse
  if((min ge 0)and(min le 9))then begin;如果分钟在0分到9分之间
    min=strjoin(['0',strtrim(min,2)])
  endif else begin
    min=strtrim(min,2)
  endelse
  ;如果秒钟在0秒到10秒之间（注意：秒钟是双精度浮点的小数）
  if((sec ge 0)and(sec lt 10))then begin
    ;把秒钟转换为无符号整型,前后的空格去掉,在前面加0
    sec=strjoin(['0',strtrim(UINT(sec),2)])
  endif else begin
    ;否则无需加0,只把秒钟转换为无符号整型,把秒钟前后的空格去掉即可
    sec=strtrim(UINT(sec),2)
  endelse
  foldername=imageprefix+year+month+day+hour+min+sec;设置文件夹名为前缀和时间
  ;新建文件夹
  file_mkdir,'C:\'+foldername+'
  Task=ENVITask('RadiometricCalibration');从ENVI任务中取得辐射校正任务
  Task.Input_Raster=RawRaster[0];设置输入的数据集：多光谱数据
  Task.Output_Data_Type='Float';设置输出的数据格式
  Task.Scale_Factor=0.1;设置比例因子
  ;设置输出地址
  Task.Output_Raster_URI='C:\'+foldername+'\'+imageprefix+'Radiance_BSQ.dat'
  Task.Execute;执行辐射校正
  DataColl=e.Data;获取数据集
  DataColl.Add,Task.Output_Raster;添加结果到数据集中
  BSQRadiance_file_address='C:\'+foldername+'\'+imageprefix+'Radiance_BSQ.dat'
  Mul_radiance_raster=e.OpenRaster(BSQRadiance_file_address)
  MulRadRaster=Mul_radiance_raster[0]
  MulRad_fid=ENVIRasterToFID(MulRadRaster)
  ;获得数据各种信息,其中增益值和偏移值由data_gains,data_offsets获得
  ENVI_FILE_QUERY,MulRad_fid,data_gains=data_gains,ns=Mul_ns,nl=Mul_nl,$
    nb=Mul_nb,dims=dims,data_offsets=data_offsets,bnames=bnames,r_fid=Mul_fid
  ;转换存储格式为BIL,指定文件夹输出,输出时在原名的基础上加了后缀'_radiance'
  Envi_doit,'convert_doit',$
    fid=MulRad_fid,$
    pos=lindgen(mul_nb),$
    dims=dims,$
    out_name='C:\'+foldername+'\'+imageprefix+'Radiance.dat',$
    o_interleave=1,$
    r_fid=BIL_fid
  ;第二部分:根据辐射亮度定标结果进行FLAASH大气校正
  radiance_file_address='C:\'+foldername+'\'+imageprefix+'Radiance.dat'
  ;通过地址得到辐射亮度文件名
  radiance_file_name=FILE_BASENAME(radiance_file_address)
  reflect_file_name=imageprefix+'Reflectance.dat';设置输出文件名
  ;输入多光谱辐射亮度数据路径,对应Input Radiance Image
  radiance_file=radiance_file_address
  ;输出大气校正结果文件路径,对应Output Relectance File
  reflect_file='C:\'+foldername+'\'+imageprefix+'Reflectance.dat'
  ;输出其他文件路径,对应 Output Directory for FLAASH Files
  qita_file='C:\'+foldername+'\
  radiance_raster=e.OpenRaster(radiance_file_address)
  ENVI_FILE_QUERY,MulRaw_fid,data_type=data_type;获取数据类型
  radiance_raster_length=Mul_RawRaster.nColumns;获取列数
  radiance_raster_width=Mul_RawRaster.nRows;获取行数
  ;获取中心波长、波长单位、波段宽度、缩放系数信息
  radiance_raster_metadata=Mul_RawRaster.Metadata
  radiance_raster_lambda=radiance_raster_metadata['WAVELENGTH'];获取中心波长
  radiance_raster_wavelength_units=radiance_raster_metadata['WAVELENGTH UNITS']
  ;获取波长单位
  radiance_raster_fwhm=radiance_raster_metadata['FWHM'];获取波段宽度
  ;获取缩放系数,此步设置为1
  radiance_raster_input_scale=MAKE_ARRAY(Mul_RawRaster.nbands,value=1,/double)
  ref=Mul_RawRaster.SPATIALREF;获取坐标信息:经纬度
  ref.ConvertFileToMap,radiance_raster_length/2,radiance_raster_width/2,MapX,MapY
  ref.ConvertMapToLonLat,MapX,MapY,radiance_raster_longitude,radiance_raster_latitude
  ;获取经纬度
  pixel_size=(ref.pixel_size)[0];获取分辨率
  ;获取时间信息,此步的时间信息是由元数据MLT.txt获取的
  ;获取时间信息字符串数组,拿参数'-T:Z'来控制分隔符
  tmpTimes=STRSPLIT(Mul_RawRaster.Time.ACQUISITION,'-T:Z',/extract)
  radiance_raster_year=FIX(tmpTimes[0])
  radiance_raster_month=FIX(tmpTimes[1])
  radiance_raster_day=FIX(tmpTimes[2])
  radiance_raster_gmt=DOUBLE(tmpTimes[3])+DOUBLE(tmpTimes[4])/60D + DOUBLE(tmpTimes[5])/60D^2
  ;通过纬度和月份来确定atmosphere_model大气模型:
  ;0-SAW;1-MLW;2-U.S. Standard;3-SAS;4-MLS;5-T
  ;第1步:根据modtran模型设置不同纬度不同月份的INT型数组
  modtran=[[0,0,0,1,1,0], $;80°~70°
    [0,0,1,1,1,0], $;70°~60°
    [1,1,1,3,3,1], $;60°~50°
    [1,1,3,3,3,3], $;50°~40°
    [3,3,3,4,4,3], $;40°~30°
    [4,4,4,5,5,4], $;30°~20°
    [5,5,5,5,5,5], $;20°~10°
    [5,5,5,5,5,5], $;10°~0°
    [5,5,5,5,5,5], $;0°~-10°
    [5,5,5,5,5,5], $;-10°~-20°
    [5,5,5,4,4,5], $;-20°~-30°
    [4,4,4,4,4,4], $;-30°~-40
    [3,3,3,3,3,3], $;-40°~-50°
    [3,3,3,1,1,3], $;-50°~-60°
    [1,1,1,1,1,1], $;-60°~-70°
    [1,1,1,1,1,1], $;-70°~-80°
    [1,1,1,0,1,1]]  ;-80°~-90°
  ;第2步:确定月份:modtran_i
  FOR i=1,6,1 DO BEGIN;月份循环
    if((radiance_raster_month eq 1)||(radiance_raster_month eq 2))then begin
      modtran_i=0
    endif else if((radiance_raster_month eq 3)||(radiance_raster_month eq 4))then begin
      modtran_i=1
    endif else if((radiance_raster_month eq 5)||(radiance_raster_month eq 6))then begin
      modtran_i=2
    endif else if((radiance_raster_month eq 7)||(radiance_raster_month eq 8))then begin
      modtran_i=3
    endif else if((radiance_raster_month eq 9)||(radiance_raster_month eq 10))then begin
      modtran_i=4
    endif else if((radiance_raster_month eq 11)||(radiance_raster_month eq 12))then begin
      modtran_i=5
    endif else begin
      a=0
    endelse
  ENDFOR
  ;第3步:确定纬度:modtran_j
  FOR j=1,17,1 DO BEGIN;纬度循环,le-小于等于,gt-大于
    if((radiance_raster_latitude le 80)&&(radiance_raster_latitude gt 70))then begin
      modtran_j=0
    endif else if((radiance_raster_latitude le 70)&&(radiance_raster_latitude gt 60))then begin
      modtran_j=1
    endif else if((radiance_raster_latitude le 60)&&(radiance_raster_latitude gt 50))then begin
      modtran_j=2
    endif else if((radiance_raster_latitude le 50)&&(radiance_raster_latitude gt 40))then begin
      modtran_j=3
    endif else if((radiance_raster_latitude le 40)&&(radiance_raster_latitude gt 30))then begin
      modtran_j=4
    endif else if((radiance_raster_latitude le 30)&&(radiance_raster_latitude gt 20))then begin
      modtran_j=5
    endif else if((radiance_raster_latitude le 20)&&(radiance_raster_latitude gt 10))then begin
      modtran_j=6
    endif else if((radiance_raster_latitude le 10)&&(radiance_raster_latitude gt 0))then begin
      modtran_j=7
    endif else if((radiance_raster_latitude le 0)&&(radiance_raster_latitude gt -10))then begin
      modtran_j=8
    endif else if((radiance_raster_latitude le -10)&&(radiance_raster_latitude gt -20))$ 
      then begin modtran_j=9
    endif else if((radiance_raster_latitude le -20)&&(radiance_raster_latitude gt -30))$
    then begin modtran_j=10
    endif else if((radiance_raster_latitude le -30)&&(radiance_raster_latitude gt -40))$
    then begin modtran_j=11
    endif else if((radiance_raster_latitude le -40)&&(radiance_raster_latitude gt -50))$
    then begin modtran_j=12
    endif else if((radiance_raster_latitude le -50)&&(radiance_raster_latitude gt -60))$
    then begin modtran_j=13
    endif else if((radiance_raster_latitude le -60)&&(radiance_raster_latitude gt -70))$
    then begin modtran_j=14
    endif else if((radiance_raster_latitude le -70)&&(radiance_raster_latitude gt -80))$
    then begin modtran_j=15
    endif else if((radiance_raster_latitude le -80)&&(radiance_raster_latitude gt -90))$
    then begin modtran_j=16
    endif else begin
      a=0
    endelse
  ENDFOR
  ;第4步:根据modtran_i,modtran_j位置确定的大气模型modtran_ij数值
  modtran_ij=modtran[modtran_i,modtran_j]
  ;获取光谱响应函数路径
  filter_func_filename='C:\Program Files\Exelis\ENVI53\resource\filterfuncs\landsat8_oli.sli'
  flaash_obj=obj_new('flaash_batch',/anc_delete);初始化FLAASH对象
  ;设置大量的输入参数
  flaash_obj->SetProperty, $
    ground_elevation=0.092, $ ;平均海拔,单位km
    ;气溶胶模型：0-No Aerosol;1-Rural;2-Maritime;3-Urban;4-Tropospheric
    aerosol_model=4, $
    ;设置FLAASH工程参数----
    hyper=0, $ ;1表示高光谱;0表示多光谱
    filter_func_filename=filter_func_filename, $设置波谱响应函数文件夹位置
    filter_func_file_index=0, $
    red_channel=4, $;0表示undefined,L8红波段为第4波段
    green_channel=3, $;0表示undefined,L8绿波段为第3波段
    blue_channel=2, $;0表示undefined,L8蓝波段为第2波段
    water_band_choice=1.13, $
    ;设置图像参数----
    data_type=data_type, $
    nspatial=radiance_raster_length, $
    nlines=radiance_raster_width, $
    margin1=0, $
    margin2=0, $
    nskip=0, $
    ;设置光谱参数
    lambda=radiance_raster_lambda, $
    wavelength_units=radiance_raster_wavelength_units, $
    fwhm=radiance_raster_fwhm, $
    input_scale=radiance_raster_input_scale, $
    ;对应FLAASH面板中的 Wavelength Recalibration,多光谱一般为0
    calc_wl_correction=0, $
    reuse_modtran_calcs=0, $
    use_square_slit_function=0, $
    convolution_method='fft', $
    ;对应 Width (number of bands) 参数,多光谱设置0即可
    polishing_res=0, $
    ;设置输入、输出文件名字及位置
    radiance_file=radiance_file, $;设置输入辐射亮度文件位置
    reflect_file=reflect_file, $;设置输出大气校正结果文件位置
    user_stem_name=imageprefix+'_', $;设置输出的其他结果文件名前缀:原先名+'_'
    modtran_directory=qita_file, $;设置输出的其他结果文件位置
    ;设置经纬度
    latitude=radiance_raster_latitude, $
    longitude=radiance_raster_longitude, $
    ;设置传感器类型、Sensor Atltitude、Ground Elevation、pixel_size
    sensor_name='Landsat-8 OLI', $
    sensor_altitude=705.0000, $ ;传感器高度
    pixel_size=pixel_size, $
    ;设置时间信息
    year=radiance_raster_year, $
    month=radiance_raster_month, $
    day=radiance_raster_day, $
    gmt=radiance_raster_gmt, $
    ;设置左下角信息
    atmosphere_model=modtran_ij, $;(根据上面获取的结果对这个进行赋值)
    water_retrieval=0, $ ;Water Retrieval参数;0表示No,1表示Yes
    water_column_multiplier=1.0000, $
    ;设置右下角信息
    aerosol_retrieval=1, $ ; 0 表示 None;1 表示 2-Band (K-T);2 表示 2-Band Over Water
    visvalue=40.0000, $ ;能见度,默认40km
    ;分别对应Multispectral Setting中Water Retrieval选项卡中的两个参数
    ;水汽反演,没有所需波段,所以设置为0,表示undefined
    water_abs_channel=0, $
    water_ref_channel=0, $
    ;对应Multispectral Setting中Kaufman-Tanre Aerosol Retrieval选项卡中的参数
    ;气溶胶反演
    kt_upper_channel=7, $ ;设置短波红外2(SWIR2)
    kt_lower_channel=4, $ ;设置红波段(Red)
    kt_cutoff=0.08, $ ;Maximum Upper Channel Reflectance
    kt_ratio=0.500, $ ;Reflectance Ratio
    cirrus_channel=0, $ ;0表示undefined
    ;对应 Advanced Setting同名参数,默认即可
    multiscatter_model=0, $
    disort_streams=8, $
    aerosol_scaleht=1.5, $;对应 Aerosol Scale Height
    co2mix=390.0000, $
    use_adjacency=1, $;对应 Use Adjacency Correction,高分辨率设置为1,低分辨率设置为0
    ;为了进行水汽反演,需要如下3个波段范围中的一个：1050-1210nm, 770-870nm, 870-1020nm
    f_resolution=15.0000, $; 而且要求此范围的波段光谱分辨率最低为15nm
    view_zenith_angle=180.0000, $
    view_azimuth=0.0000, $
    use_tiling=0, $ ;1-Yes;0-No
    tile_size=100.0000, $
    output_scale=1;输出结果缩放系数,输出结果放大了1倍,变为UINT数据类型
  radiance_raster.Close;执行FLAASH之前,必须在ENVI中把输入文件关闭
  flaash_obj->processImage;开始执行FLAASH
  ;第三部分:大气辐射传输方程法的波段运算
  RefelectanceRasters=e.OpenRaster(reflect_file)
  RefelectanceRaster=RefelectanceRasters[0]
  fid=ENVIRasterToFID(RefelectanceRaster)
  Mul_map_info=envi_get_map_info(fid=fid)
  ENVI_FILE_QUERY,fid,DIMS=dims
  ;计算NDVI归一化植被指数
  envi_doit, 'math_doit', $
    FID = [fid, fid], $
    pos=[3,4], $
    dims=dims, $
    exp='(float(b5)-float(b4)) / (float(b5)+float(b4)+0.000001)',$
    out_name='C:\'+foldername+'\'+imageprefix+'NDVI.dat',$
    r_fid=NDVI_fid
  ndvi=ENVI_GET_DATA(FID=NDVI_fid,DIMS=DIMS,POS=0)
  ;计算Fv(或者说Pv):植被在混合像元内的构成比例即植被覆盖度;Fv=(NDVI-NDVIs)/(NDVIv-NDVIs)
  ;NDVIsoil为完全是裸土或无植被覆盖区域的NDVI值
  ;NDVIveg则代表完全被植被所覆盖的像元的NDVI值,即纯植被像元的NDVI值
  ;取经验值NDVIv=0.70和NDVIs=0.05
  ;即当某个像元的NDVI大于0.70时,Fv取值为1;当NDVI小于0.05,Fv取值为0
  ;这个经验值应该是Sobrino提出的
  FV=((NDVI GT 0.7)*1+(NDVI LT 0.05)*0+(NDVI GE 0.05 AND NDVI LE 0.7)*((NDVI-0.05)/(0.7-0.05)))^2
  ;输出Fv到硬盘
  envi_write_envi_file,Fv,$
    out_name='C:\'+foldername+'\'+imageprefix+'Fv.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Fv_fid,$
    map_info=Mul_map_info
  ;计算自然表面在L8B10波长区间内的辐射率(ENVI-IDL的技术殿堂,ETM+数据)
  emission_surface=0.9625 + 0.0614*fv - 0.0461*fv^2
  envi_write_envi_file,emission_surface,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission_surface.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=es_fid,$
    map_info=Mul_map_info
  ;计算城镇像元在L8B10波长区间内的辐射率(ENVI-IDL的技术殿堂,ETM+数据)
  emission_building=0.9589 + 0.086*fv - 0.0671*fv^2
  envi_write_envi_file,emission_building,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission_building.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=eb_fid,$
    map_info=Mul_map_info
  ;计算地表发射率(ENVI-IDL的技术殿堂,ETM+数据)
  ;将遥感影像分为水体、城镇和自然表面3种类型;水体像元的比辐射率赋值为0.995
  emission=(ndvi le 0)*0.995+(ndvi gt 0 and ndvi lt 0.7)*emission_building+(ndvi ge 0.7)*emission_surface
  envi_write_envi_file,emission,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=e_fid,$
    map_info=Mul_map_info
  ;对热红外波段进行辐射定标
  Task = ENVITask('RadiometricCalibration')
  Task.Input_Raster = RawRaster[3] ;热红外波段
  Task.Output_Data_Type = 'Float'
  Task.Scale_Factor = 1 ;设置比例因子
  Task.Output_Raster_URI = 'C:\'+foldername+'\'+imageprefix+'Radiance_Tir.dat'
  Task.Execute
  DataColl = e.Data
  DataColl.Add, Task.Output_Raster
  ;取得辐射定标结果,求解B(Ts)
  TirRadiance_file_address='C:\'+foldername+'\'+imageprefix+'Radiance_Tir.dat'
  Tir_radiance_raster=e.OpenRaster(TirRadiance_file_address)
  TirRad_fid=ENVIRasterToFID(Tir_radiance_raster)
  Radiance_B10=envi_get_data(fid=TirRad_fid,dims=dims,pos=0)
  upwelling_radiance=0.28;上行辐射
  downwelling_radiance=0.53;下行辐射
  transmission=0.96;大气透过率
  ev=0.98672;单窗法中使用,植物在L8B10波长区间内的辐射率
  em=0.96767;单窗法中使用,建筑表面的辐射率
  a6=-62.7182;单窗法中使用,系数,0-70摄氏度下为-67.355351
  b6=0.4339;单窗法中使用,系数,0-70摄氏度下为0.458606
  T=12;单窗法中使用,卫星过境时的近地表温度(摄氏度)
  T0=273+T;单窗法中使用,卫星过境时的近地表温度(开尔文)
  Ta=17.9769+0.91715*T0;单窗法中使用,大气平均作用温度
  ev10=0.98672
  ev11=0.98990
  es10=0.96767
  es11=0.97790
  a10=-64.4661
  a11=-68.8678
  b10=0.4398
  b11=0.4755
  ;B(Ts)是用Planck函数表示的黑体热辐射强度,或者说温度为T的黑体在热红外波段的辐射亮度值
  ;这个公式是由大气辐射传输方程转化而来的;L=Lu+B(Ts)*e*t+(1-e)*t*Ld
  BTs=(Radiance_B10-upwelling_radiance-transmission*(1-emission)*downwelling_radiance) / (transmission*emission)
  envi_write_envi_file,BTs,$
    out_name='C:\'+foldername+'\'+imageprefix+'B(Ts).dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=BTs_fid,$
    map_info=Mul_map_info
  ;输出地表温度;这里输出的是摄氏度
  ;这是Planck函数的近似式;T=K2/(alog(K1/B(Ts)+1))
  ;K1=h*c/(k*lamda^5)=774.89,k2=2*h*c^2/(lamda^5)=1321.08
  Tdq=1321.08/(alog(774.89/BTs+1))-273
  envi_write_envi_file,Tdq,$
    out_name='C:\'+foldername+'\'+imageprefix+'Tdq.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Tdq_fid,$
    map_info=Mul_map_info
  ;第四部分:单窗法的波段运算
  ;计算Rv植被的温度比率(覃志豪,TM6)
  Rv=0.9332+0.0585*fv
  envi_write_envi_file,Rv,$
    out_name='C:\'+foldername+'\'+imageprefix+'Rv.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Rv_fid,$
    map_info=Mul_map_info
  ;计算Rm建筑表面的温度比率(覃志豪,TM6)
  Radiom=0.9886+0.1287*fv
  envi_write_envi_file,Radiom,$
    out_name='C:\'+foldername+'\'+imageprefix+'Rm.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Rm_fid,$
    map_info=Mul_map_info
  ;计算de;delta_emission,发射率的微分(覃志豪,TM6)
  de=(fv le 0.5)*0.003796*fv+(fv gt 0.5)*0.003796*(1-fv)
  envi_write_envi_file,de,$
    out_name='C:\'+foldername+'\'+imageprefix+'de.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=de_fid,$
    map_info=Mul_map_info
  ;做Landsat8 热红外波段的亮温
  Task = ENVITask('RadiometricCalibration')
  Task.CALIBRATION_TYPE = 'Brightness Temperature';校正类型为亮度温度
  Task.Input_Raster = RawRaster[3] ;第三个数据集是热红外波段
  Task.Output_Data_Type = 'Float'
  Task.Scale_Factor = 1
  Task.Output_Raster_URI = 'C:\'+foldername+'\'+imageprefix+'BrightnessTemperature.dat'
  Task.Execute
  DataColl = e.Data
  DataColl.Add, Task.Output_Raster
  ;取得亮温结果
  BrightnessTemperature_file_address='C:\'+foldername+'\'+imageprefix+'BrightnessTemperature.dat'
  Tir_BT_raster=e.OpenRaster(BrightnessTemperature_file_address)
  TirBT_fid=ENVIRasterToFID(Tir_BT_raster)
  T6=envi_get_data(fid=TirBT_fid,dims=dims,pos=0)
  T11=envi_get_data(fid=TirBT_fid,dims=dims,pos=1)
  ;epsilon是地表发射率(覃志豪,TM6)
  epsilon=fv*rv*ev+(1-fv)*radiom*em+de
  envi_write_envi_file,epsilon,$
    out_name='C:\'+foldername+'\'+imageprefix+'epsilon.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=epsilon_fid,$
    map_info=Mul_map_info
  ;计算统计信息
  stats = ENVIRasterStatistics(Tir_BT_raster,/COVARIANCE)
  ;协方差矩阵的第二个值是波段10和波段11的协方差
  COVARIANCE=stats.COVARIANCE[1]
  ;协方差矩阵的第一个值是波段10方差
  VARIANCE_B10=stats.COVARIANCE[0]
  ;下面是求大气水汽的公式
  wv=-9.674*(COVARIANCE/VARIANCE_B10)*(COVARIANCE/VARIANCE_B10)+0.653*(COVARIANCE/VARIANCE_B10)+9.087
  CASE 1 OF
    (wv ge 0.2)and(wv le 2.0):transmission10=0.9220-0.0780*wv
    (wv gt 2.0)and(wv le 5.6):transmission10=1.0222-0.1330*wv
    (wv ge 5.6)and(wv le 6.8):transmission10=0.5422-0.0440*wv
    (wv gt 6.8):print,'大气水汽不在合理之内'
    (wv lt 0.2):print,'大气水汽不在合理之内'
  ENDCASE
  ;C6是为了简化公式用的中间变量
  C6=epsilon*transmission10
  envi_write_envi_file,C6,$
    out_name='C:\'+foldername+'\'+imageprefix+'C6.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=C6_fid,$
    map_info=Mul_map_info
  ;D6是为了简化公式用的中间变量
  D6=(1-transmission10)*(1+(1-epsilon)*transmission10)
  envi_write_envi_file,D6,$
    out_name='C:\'+foldername+'\'+imageprefix+'D6.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=D6_fid,$
    map_info=Mul_map_info
  ;这就是单窗法(覃志豪,TM6)的公式,这里输出的是摄氏度
  Tmw=(a6*(1-C6-D6)+(b6*(1-C6-D6)+C6+D6)*T6-D6*Ta)/C6-273
  envi_write_envi_file,Tmw,$
    out_name='C:\'+foldername+'\'+imageprefix+'Tmw.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Tsurface_fid,$
    map_info=Mul_map_info
  Rs=0.9902+0.1068*fv
  envi_write_envi_file,Rs,$
    out_name='C:\'+foldername+'\'+imageprefix+'Rs.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Rs_fid,$
    map_info=Mul_map_info
  e10=fv*rv*ev10+(1-fv)*rs*es10+de
  envi_write_envi_file,e10,$
    out_name='C:\'+foldername+'\'+imageprefix+'e10.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=e10_fid,$
    map_info=Mul_map_info
  e11=fv*rv*ev11+(1-fv)*rs*es11+de
  envi_write_envi_file,e11,$
    out_name='C:\'+foldername+'\'+imageprefix+'e11.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=e11_fid,$
    map_info=Mul_map_info
  tau10=-0.1146*wv+1.0286
  tau11=-0.1568*wv+1.0083
  C10=e10*tau10
  envi_write_envi_file,C10,$
    out_name='C:\'+foldername+'\'+imageprefix+'C10.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=C10_fid,$
    map_info=Mul_map_info
  C11=e11*tau11
  envi_write_envi_file,C11,$
    out_name='C:\'+foldername+'\'+imageprefix+'C11.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=C11_fid,$
    map_info=Mul_map_info
  D10=(1-tau10)*(1+(1-e10)*tau10)
  envi_write_envi_file,D10,$
    out_name='C:\'+foldername+'\'+imageprefix+'D10.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=D10_fid,$
    map_info=Mul_map_info
  D11=(1-tau11)*(1+(1-e11)*tau11)
  envi_write_envi_file,D11,$
    out_name='C:\'+foldername+'\'+imageprefix+'D11.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=D11_fid,$
    map_info=Mul_map_info
  E0=D11*C10-D10*C11
  envi_write_envi_file,E0,$
    out_name='C:\'+foldername+'\'+imageprefix+'E0.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=E0_fid,$
    map_info=Mul_map_info
  E1=D11*(1-C10-D10)/E0
  envi_write_envi_file,E1,$
    out_name='C:\'+foldername+'\'+imageprefix+'E1.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=E1_fid,$
    map_info=Mul_map_info
  E2=D10*(1-C11-D11)/E0
  envi_write_envi_file,E2,$
    out_name='C:\'+foldername+'\'+imageprefix+'E2.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=E2_fid,$
    map_info=Mul_map_info
  A0=E1*a10+E2*a11
  envi_write_envi_file,A0,$
    out_name='C:\'+foldername+'\'+imageprefix+'A0.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=A0_fid,$
    map_info=Mul_map_info
  A=D10/E0
  envi_write_envi_file,A,$
    out_name='C:\'+foldername+'\'+imageprefix+'A.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=A_fid,$
    map_info=Mul_map_info
  A1=1+A+E1*b10
  envi_write_envi_file,A1,$
    out_name='C:\'+foldername+'\'+imageprefix+'A1.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=A1_fid,$
    map_info=Mul_map_info
  A2=A+E2*b11
  envi_write_envi_file,A2,$
    out_name='C:\'+foldername+'\'+imageprefix+'A2.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=A2_fid,$
    map_info=Mul_map_info
  Tsw=A0+A1*T6-A2*T11-273
  envi_write_envi_file,Tsw,$
    out_name='C:\'+foldername+'\'+imageprefix+'Tsw.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Tsurface_fid,$
    map_info=Mul_map_info
END