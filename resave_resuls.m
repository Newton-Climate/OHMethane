
resave('case0.mat')
resave('case1.mat')
resave('case2.mat')
resave('case3.mat')
  resave('case4.mat')
  resave('case5.mat')


function ems_struct = ems2struct(ems)
  ems_struct.nh_ch4    = ems(:,1);
  ems_struct.sh_ch4    = ems(:,2);
  ems_struct.nh_ch4c13 = ems(:,3);
  ems_struct.sh_ch4c13 = ems(:,4);
  ems_struct.nh_mcf    = ems(:,5);
  ems_struct.sh_mcf    = ems(:,6);
  ems_struct.nh_n2o    = ems(:,7);
  ems_struct.sh_n2o    = ems(:,8);
  ems_struct.nh_c2h6   = ems(:,9);
  ems_struct.sh_c2h6   = ems(:,10);
  ems_struct.nh_oh     = ems(:,11);
  ems_struct.sh_oh     = ems(:,12);
  ems_struct.nh_co     = ems(:,13);
  ems_struct.sh_co     = ems(:,14);
  ems_struct.tau_TS    = ems(:,15);
  ems_struct.kX_NH     = ems(:,16);
ems_struct.kX_SH     = ems(:,17);
end


function out = resave(filename)
  results = open(filename)
  ems = ems2struct(results.ems_anal);
  con = results.out_anal;
  t = results.St;
save(filename);
out = true;
end
