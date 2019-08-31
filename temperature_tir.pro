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
