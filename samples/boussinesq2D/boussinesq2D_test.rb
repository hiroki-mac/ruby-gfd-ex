# 2次元ブシネスク方程式
require "numru/gphys"
require "numru/ganalysis"
require "../../lib/timeintegration.rb"
require "../../lib/gphysmath.rb"
require "numru/ganalysis/ganalysis_ext"
require "pry"
include NumRu
include NMath
include GPMath

# 定数設定
dt = UNumeric[ARGV[0].to_f,"s"]   # 時間刻み幅
HD = UNumeric[5E4.to_f,"m2.s-1"]  # 拡散係数 (∇^2)
T = UNumeric[100,"s"]         # 積分時間

TempT = UNumeric[300,"K"]         # 上部境界温度
TempB = UNumeric[305,"K"]         # 下部境界温度
Temp0 = (TempT + TempB)/2         # 平均温度
Depth = UNumeric[1000,"m"]        # 領域の厚さ
Grav  = UNumeric[9.8,"m.s-2"]     # 重力加速度
Nu    = UNumeric[1.0E2,"m2.s-1"]  # 動粘性係数
Kappa = UNumeric[1.0E2,"m2.s-1"]  # 熱拡散係数

# 無次元数
Ra = ( Grav*(TempB - TempT)*Depth**3 )/(Kappa * Nu * Temp0)
Pr = Nu/Kappa
print "Ra = #{Ra}\n"
print "Pr = #{Pr}\n"

# 領域設定
# 2次元平面 1 km x 1 km. x方向は周期境界. y方向は非周期境界.
grid = GPhys.domain2D(
 {"range"=>0..8E3, "grid_points"=>256, "cyclic"=> true, "units"=>"m", "virtual_points"=>0},
 {"range"=>0..1E3, "grid_points"=>33, "cyclic"=> false, "units"=>"m", "cell_center"=>false} )
 
# 変数設定
# 渦度
vor = GPhys.make_var(grid,
  {"name"=>"vor", "long_name"=>"vorticity", "units"=>"s-1"} )
# 流線関数
stf = GPhys.make_var(grid,
  {"name"=>"stf", "long_name"=>"streamfunction", "units"=>"s-1"} )
# x方向速度
u = GPhys.make_var(grid,
  {"name"=>"u", "long_name"=>"u-velocity", "units"=>"m.s-1"} )
# y方向速度
v = GPhys.make_var(grid,
  {"name"=>"v", "long_name"=>"v-velocity", "units"=>"m.s-1"} )
# 温度
temp = GPhys.make_var(grid,
  {"name"=>"Temp", "long_name"=>"Temperature", "units"=>"K"} )
# 密度
rho = GPhys.make_var(grid,
  {"name"=>"rho", "long_name"=>"density perturbation", "units"=>"kg.m-3"} )


##############
# 初期分布設定 #
##############
#temp = temp*0 + Temp0
#temp[true, 0, false] = temp[true, 0, false]*0 + TempB
#temp[true,-1, false] = temp[true,-1, false]*0 + TempT
#temp[31..33, 0, false] = temp[31..33, 0, false] + 1.0
jm = grid.shape[1]
jm.times{|j|
  temp[true,j,false] = temp[true,j,false]*0 + (TempB*(jm-1-j) + TempT*(j) )/(jm-1)
}
temp_base = temp.copy
temp = temp*0.0
#temp[32,1..32,false] = temp[32,0..31,false]*1.00
temp[32,1..32,false] = temp[32,0..31,false] + 0.1

# test data
# im = grid.shape[0]
# jm.times{|j|
#   im.times{|i|
#   temp[i,j,false] = temp[i,j,false]*0 + sin(2*PI * i.to_f/(im)) * sin(1*PI * j.to_f/(jm))
#   }
# }


# # グリッドデータで与える場合
# xc0 = 4000E3; yc0 = 5000E3; xc1 = 6000E3; yc1 = 5000E3; rlen = 1000E3
# # xc0 = 150; yc0 = 0; xc1 = 210; yc1 = 0; rlen = 20
# xval = vor.coord(0).val; yval = vor.coord(1).val
# val = vor.val
# grid.shape[1].times{|j|
#   grid.shape[0].times{|i|
#     r0 = sqrt( (xval[i]-xc0)**2 + (yval[j]-yc0)**2 )
#     r1 = sqrt( (xval[i]-xc1)**2 + (yval[j]-yc1)**2 )
#     val[i,j,0] = (exp((rlen-r0)/rlen*2) + exp((rlen-r1)/rlen*2))*1.0E-5
#   }
# }
# vor.replace_val(val)

GPhys.prepare_fft_wavenumber(vor)

#binding.pry

# Vorticity equation on double cyclic 2D-plane
boussinesq2D = lambda {|gpna,t|
  GPhys.use_fft
  # 変数格納（渦度、温度）
  zeta = gpna[0]; temp = gpna[1]
  # 流線関数の計算 ζ = △(ψ) → ψ = △^-1(ζ)
  psi = invlap(zeta)
  #psi = psi-psi[true,0,false] # 境界条件
  zeta_tendency = -jacobian(psi,zeta) + (1.0/Temp0)*Grav*temp.dx + Nu*lap(zeta)
  temp_tendency = -jacobian(psi,temp) - psi.dx*(TempT-TempB)/Depth + Kappa*lap(temp)
  return NArray.to_na([zeta_tendency, temp_tendency])
}

ofn_vor = "vor.nc"
ofn_stf = "stf.nc"
#ofn_u = "u.nc"
#ofn_v = "v.nc"
ofn_temp = "Temp.nc"
#ofn_rho = "rho.nc"

gpna = NArray.to_na([vor, temp])

i = 0
# 時間発展ループ
while (dt*i < T) do
  if (i%20 == 0) then
    auto_write(ofn_vor, gpna[0], i, false)
    auto_write(ofn_temp, gpna[1]+temp_base, i, false)
    stf.replace_val(invlap(gpna[0]).val)
    auto_write(ofn_stf, stf, i, false)
  end
    
  gpna = TimeIntegration.rk3(gpna,dt,1,&boussinesq2D)
  # 境界条件
  gpna[0][true, 0,false] = gpna[0][true, 0,false]*0
  gpna[0][true,-1,false] = gpna[0][true,-1,false]*0
  gpna[1][true, 0,false] = gpna[1][true, 0,false]*0
  gpna[1][true,-1,false] = gpna[1][true,-1,false]*0
#  gpna[1][31..33,   0,false] = gpna[1][31..33,   0,false] + 1.0

  p i*dt if (i%20 == 0)
  i += 1
  if gpna[0].max.to_f.nan? then
    raise "!!! Numerical Instability occurs !!!"
  end
end
