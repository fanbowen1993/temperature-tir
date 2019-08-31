;��������:�������䴫�䷽�̷��Ĳ�������
  RefelectanceRasters=e.OpenRaster(reflect_file)
  RefelectanceRaster=RefelectanceRasters[0]
  fid=ENVIRasterToFID(RefelectanceRaster)
  Mul_map_info=envi_get_map_info(fid=fid)
  ENVI_FILE_QUERY,fid,DIMS=dims
  ;����NDVI��һ��ֲ��ָ��
  envi_doit, 'math_doit', $
    FID = [fid, fid], $
    pos=[3,4], $
    dims=dims, $
    exp='(float(b5)-float(b4)) / (float(b5)+float(b4)+0.000001)',$
    out_name='C:\'+foldername+'\'+imageprefix+'NDVI.dat',$
    r_fid=NDVI_fid
  ndvi=ENVI_GET_DATA(FID=NDVI_fid,DIMS=DIMS,POS=0)
  ;����Fv(����˵Pv):ֲ���ڻ����Ԫ�ڵĹ��ɱ�����ֲ�����Ƕ�;Fv=(NDVI-NDVIs)/(NDVIv-NDVIs)
  ;NDVIsoilΪ��ȫ����������ֲ�����������NDVIֵ
  ;NDVIveg�������ȫ��ֲ�������ǵ���Ԫ��NDVIֵ,����ֲ����Ԫ��NDVIֵ
  ;ȡ����ֵNDVIv=0.70��NDVIs=0.05
  ;����ĳ����Ԫ��NDVI����0.70ʱ,FvȡֵΪ1;��NDVIС��0.05,FvȡֵΪ0
  ;�������ֵӦ����Sobrino�����
  FV=((NDVI GT 0.7)*1+(NDVI LT 0.05)*0+(NDVI GE 0.05 AND NDVI LE 0.7)*((NDVI-0.05)/(0.7-0.05)))^2
  ;���Fv��Ӳ��
  envi_write_envi_file,Fv,$
    out_name='C:\'+foldername+'\'+imageprefix+'Fv.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Fv_fid,$
    map_info=Mul_map_info
  ;������Ȼ������L8B10���������ڵķ�����(ENVI-IDL�ļ�������,ETM+����)
  emission_surface=0.9625 + 0.0614*fv - 0.0461*fv^2
  envi_write_envi_file,emission_surface,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission_surface.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=es_fid,$
    map_info=Mul_map_info
  ;���������Ԫ��L8B10���������ڵķ�����(ENVI-IDL�ļ�������,ETM+����)
  emission_building=0.9589 + 0.086*fv - 0.0671*fv^2
  envi_write_envi_file,emission_building,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission_building.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=eb_fid,$
    map_info=Mul_map_info
  ;����ر�����(ENVI-IDL�ļ�������,ETM+����)
  ;��ң��Ӱ���Ϊˮ�塢�������Ȼ����3������;ˮ����Ԫ�ıȷ����ʸ�ֵΪ0.995
  emission=(ndvi le 0)*0.995+(ndvi gt 0 and ndvi lt 0.7)*emission_building+(ndvi ge 0.7)*emission_surface
  envi_write_envi_file,emission,$
    out_name='C:\'+foldername+'\'+imageprefix+'emission.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=e_fid,$
    map_info=Mul_map_info
  ;���Ⱥ��Ⲩ�ν��з��䶨��
  Task = ENVITask('RadiometricCalibration')
  Task.Input_Raster = RawRaster[3] ;�Ⱥ��Ⲩ��
  Task.Output_Data_Type = 'Float'
  Task.Scale_Factor = 1 ;���ñ�������
  Task.Output_Raster_URI = 'C:\'+foldername+'\'+imageprefix+'Radiance_Tir.dat'
  Task.Execute
  DataColl = e.Data
  DataColl.Add, Task.Output_Raster
  ;ȡ�÷��䶨����,���B(Ts)
  TirRadiance_file_address='C:\'+foldername+'\'+imageprefix+'Radiance_Tir.dat'
  Tir_radiance_raster=e.OpenRaster(TirRadiance_file_address)
  TirRad_fid=ENVIRasterToFID(Tir_radiance_raster)
  Radiance_B10=envi_get_data(fid=TirRad_fid,dims=dims,pos=0)
  upwelling_radiance=0.28;���з���
  downwelling_radiance=0.53;���з���
  transmission=0.96;����͸����
  ev=0.98672;��������ʹ��,ֲ����L8B10���������ڵķ�����
  em=0.96767;��������ʹ��,��������ķ�����
  a6=-62.7182;��������ʹ��,ϵ��,0-70���϶���Ϊ-67.355351
  b6=0.4339;��������ʹ��,ϵ��,0-70���϶���Ϊ0.458606
  T=12;��������ʹ��,���ǹ���ʱ�Ľ��ر��¶�(���϶�)
  T0=273+T;��������ʹ��,���ǹ���ʱ�Ľ��ر��¶�(������)
  Ta=17.9769+0.91715*T0;��������ʹ��,����ƽ�������¶�
  ev10=0.98672
  ev11=0.98990
  es10=0.96767
  es11=0.97790
  a10=-64.4661
  a11=-68.8678
  b10=0.4398
  b11=0.4755
  ;B(Ts)����Planck������ʾ�ĺ����ȷ���ǿ��,����˵�¶�ΪT�ĺ������Ⱥ��Ⲩ�εķ�������ֵ
  ;�����ʽ���ɴ������䴫�䷽��ת��������;L=Lu+B(Ts)*e*t+(1-e)*t*Ld
  BTs=(Radiance_B10-upwelling_radiance-transmission*(1-emission)*downwelling_radiance) / (transmission*emission)
  envi_write_envi_file,BTs,$
    out_name='C:\'+foldername+'\'+imageprefix+'B(Ts).dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=BTs_fid,$
    map_info=Mul_map_info
  ;����ر��¶�;��������������϶�
  ;����Planck�����Ľ���ʽ;T=K2/(alog(K1/B(Ts)+1))
  ;K1=h*c/(k*lamda^5)=774.89,k2=2*h*c^2/(lamda^5)=1321.08
  Tdq=1321.08/(alog(774.89/BTs+1))-273
  envi_write_envi_file,Tdq,$
    out_name='C:\'+foldername+'\'+imageprefix+'Tdq.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Tdq_fid,$
    map_info=Mul_map_info
  ;���Ĳ���:�������Ĳ�������
  ;����Rvֲ�����¶ȱ���(��־��,TM6)
  Rv=0.9332+0.0585*fv
  envi_write_envi_file,Rv,$
    out_name='C:\'+foldername+'\'+imageprefix+'Rv.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Rv_fid,$
    map_info=Mul_map_info
  ;����Rm����������¶ȱ���(��־��,TM6)
  Radiom=0.9886+0.1287*fv
  envi_write_envi_file,Radiom,$
    out_name='C:\'+foldername+'\'+imageprefix+'Rm.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=Rm_fid,$
    map_info=Mul_map_info
  ;����de;delta_emission,�����ʵ�΢��(��־��,TM6)
  de=(fv le 0.5)*0.003796*fv+(fv gt 0.5)*0.003796*(1-fv)
  envi_write_envi_file,de,$
    out_name='C:\'+foldername+'\'+imageprefix+'de.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=de_fid,$
    map_info=Mul_map_info
  ;��Landsat8 �Ⱥ��Ⲩ�ε�����
  Task = ENVITask('RadiometricCalibration')
  Task.CALIBRATION_TYPE = 'Brightness Temperature';У������Ϊ�����¶�
  Task.Input_Raster = RawRaster[3] ;���������ݼ����Ⱥ��Ⲩ��
  Task.Output_Data_Type = 'Float'
  Task.Scale_Factor = 1
  Task.Output_Raster_URI = 'C:\'+foldername+'\'+imageprefix+'BrightnessTemperature.dat'
  Task.Execute
  DataColl = e.Data
  DataColl.Add, Task.Output_Raster
  ;ȡ�����½��
  BrightnessTemperature_file_address='C:\'+foldername+'\'+imageprefix+'BrightnessTemperature.dat'
  Tir_BT_raster=e.OpenRaster(BrightnessTemperature_file_address)
  TirBT_fid=ENVIRasterToFID(Tir_BT_raster)
  T6=envi_get_data(fid=TirBT_fid,dims=dims,pos=0)
  T11=envi_get_data(fid=TirBT_fid,dims=dims,pos=1)
  ;epsilon�ǵر�����(��־��,TM6)
  epsilon=fv*rv*ev+(1-fv)*radiom*em+de
  envi_write_envi_file,epsilon,$
    out_name='C:\'+foldername+'\'+imageprefix+'epsilon.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=epsilon_fid,$
    map_info=Mul_map_info
  ;����ͳ����Ϣ
  stats = ENVIRasterStatistics(Tir_BT_raster,/COVARIANCE)
  ;Э�������ĵڶ���ֵ�ǲ���10�Ͳ���11��Э����
  COVARIANCE=stats.COVARIANCE[1]
  ;Э�������ĵ�һ��ֵ�ǲ���10����
  VARIANCE_B10=stats.COVARIANCE[0]
  ;�����������ˮ���Ĺ�ʽ
  wv=-9.674*(COVARIANCE/VARIANCE_B10)*(COVARIANCE/VARIANCE_B10)+0.653*(COVARIANCE/VARIANCE_B10)+9.087
  CASE 1 OF
    (wv ge 0.2)and(wv le 2.0):transmission10=0.9220-0.0780*wv
    (wv gt 2.0)and(wv le 5.6):transmission10=1.0222-0.1330*wv
    (wv ge 5.6)and(wv le 6.8):transmission10=0.5422-0.0440*wv
    (wv gt 6.8):print,'����ˮ�����ں���֮��'
    (wv lt 0.2):print,'����ˮ�����ں���֮��'
  ENDCASE
  ;C6��Ϊ�˼򻯹�ʽ�õ��м����
  C6=epsilon*transmission10
  envi_write_envi_file,C6,$
    out_name='C:\'+foldername+'\'+imageprefix+'C6.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=C6_fid,$
    map_info=Mul_map_info
  ;D6��Ϊ�˼򻯹�ʽ�õ��м����
  D6=(1-transmission10)*(1+(1-epsilon)*transmission10)
  envi_write_envi_file,D6,$
    out_name='C:\'+foldername+'\'+imageprefix+'D6.dat',$
    o_interleave=1,$
    offset=Mul_offset,$
    r_fid=D6_fid,$
    map_info=Mul_map_info
  ;����ǵ�����(��־��,TM6)�Ĺ�ʽ,��������������϶�
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
