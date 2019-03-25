require "numru/gphys"
require "numru/ganalysis"
require "../../lib/timeintegration.rb"
require "../../lib/gphysmath.rb"

require "pry"
include NumRu
include NMath
include GPMath
include SphericalHarmonics
Pry.config.pager = false
#require "profile"

dt = UNumeric[ARGV[0].to_f,"s"] #0.01 # 時間刻み幅
p dt
HD = UNumeric[1E5.to_f,"m2.s-1"]  # 拡散係数 (∇^2)
p HD

C = UNumeric[1.0E-1,"m.s-1"] # 移流速度

T = UNumeric[86400*1000,"s"] # 積分時間

Beta = UNumeric[0.0E-11, "m-1.s-1" ]

R = GAnalysis::Planet.radius
Omega = GAnalysis::Planet.omega*0
g = GAnalysis::Met.g*0
D2R = PI/180.0


# 領域設定
# 周期境界の2次元平面
#grid = GPhys.domain2D(
#  {"range"=>0..1E7, "grid_points"=>64, "cyclic"=> true, "units"=>"m"},
#  {"range"=>0..1E7, "grid_points"=>64, "cyclic"=> true, "units"=>"m"} )

# 2次元球面  
grid = GPhys.domain2Dsphere(64,R) # 引数は経度方向の格子点数
  
# 変数設定（渦度、発散、表面変位）
vor = GPhys.make_var(grid,
  {"name"=>"vor", "long_name"=>"vorticity", "units"=>"s-1"} )
div = GPhys.make_var(grid,
  {"name"=>"div", "long_name"=>"divergence", "units"=>"s-1"} )
h = GPhys.make_var(grid,
  {"name"=>"h", "long_name"=>"height-disturbance", "units"=>"m"} )

dx = vor.coord(0).val[1]-vor.coord(0).val[0]
p dx

# 球面調和関数変換の準備
GPhys.prepare_spht(vor)
GPhys.use_spht
coslat2 = SphericalHarmonics.coslat()*SphericalHarmonics.coslat()
coslat = SphericalHarmonics.coslat()
f = SphericalHarmonics.sinlat()*2.0*Omega

##############
# 初期分布設定 #
##############
# グリッドデータで与える場合
xval = vor.coord(0).val; yval = vor.coord(1).val
val = vor.val
#xc0 = PI-PI/4; yc0 = PI; xc1 = PI+PI/4; yc1 = PI; rlen = PI/8
#xc0 = 4000E3; yc0 = 5000E3; xc1 = 6000E3; yc1 = 5000E3; rlen = 1000E3
xc0 = 150; yc0 = 0; xc1 = 210; yc1 = 0; rlen = 20
grid.shape[1].times{|j|
  grid.shape[0].times{|i|
    r0 = sqrt( (xval[i]-xc0)**2 + (yval[j]-yc0)**2 )
    r1 = sqrt( (xval[i]-xc1)**2 + (yval[j]-yc1)**2 )
#    val[i,j,0] = (exp((rlen-r0)/rlen*2) + exp((rlen-r1)/rlen*2))*1.0E-5
    val[i,j,0] = sin(xval[i]*D2R)*cos(yval[j]*D2R)*1.0E-5 # Y^1_1
  }
}
vor.replace_val(val)

# スペクトルで与える場合
# vor_sp = vor.sp 
# vor_sp(m,n) は Y^m_n の係数, m: 東西波数, n:全波数
# binding.pry
# vor_sp[0..5,5..7,false] = NArray.dcomplex(6,3,1).randomn*1E-6
# vor = vor_sp.gp

# Williamson et al. (1992) テスト1 極を越えるcos型の山の移流
val = vor.val; h_val = h.val
alpha = PI/2; yc = 0*D2R;
grid.shape[1].times{|j|
  grid.shape[0].times{|i|
    val[i,j,0] = -(sin(yval[j]*D2R)*cos(alpha) - cos(xval[i]*D2R)*cos(yval[j]*D2R)*sin(alpha))*R*(-50.0)
    r = acos(sin(yval[j]*D2R)*sin(yc) + cos(yval[j]*D2R)*cos(xval[i]*D2R-1.5*PI)*cos(yc))*R
    if (r < R/3) then
      h_val[i,j,0] = (1000.0/2.0)*((PI*r/(R/3)).cos + UNumeric[1.0,"1"])
    else
      h_val[i,j,0] = 0.0 #
    end   
  }
}
vor.replace_val(val); h.replace_val(h_val)

vor = lap(vor).set_att("units","s-1").set_att("long_name","vorticity").rename("vor")


# Shallow water equation on a double cyclic 2D-plane
# shallow_2Dplane = lambda {|gpna,t|
#   GPhys.use_fft
#   zeta = gpna[0]
#   bc = 2 # cyclic boundary condition
#   psi = invlap(zeta)
#   zeta_tendency = -jacobian(psi,zeta) - Beta*psi.dx + HD*lap(zeta)
#   return NArray.to_na([zeta_tendency])
# }

# Shallow water equation on a sphere
shallow_sphere = lambda {|gpna,t|
  GPhys.output_in_spectral
  vor = gpna[0] # 渦度
  div = gpna[1] # 発散
  h = gpna[2].gp   # 高さ 
  psi = invlap(vor) # 流線関数
  chi = invlap(div) # 速度ポテンシャル
  u_cos = -psi.cos2_dy.gp + chi.dx.gp # u*cos(lat)
  v_cos = chi.cos2_dy.gp + psi.dx.gp  # v*cos(lat)    

  en = (u_cos*u_cos + v_cos*v_cos)/(2.0*coslat2) # 運動エネルギー
  abs_vor = vor.gp + f # 絶対渦度

  # 時間変化率
  vor_dt = ( -(abs_vor*u_cos/coslat2).dx - (abs_vor*v_cos).dy + HD*(lap(vor) + 2/(R*R)*vor) )
  div_dt = (  (abs_vor*v_cos/coslat2).dx - (abs_vor*u_cos).dy - (lap(h*g + en)) + HD*(lap(div) + 2/(R*R)*div) )
  h_dt   = ( -(h*u_cos/coslat2).dx - (v_cos*h).dy + HD*lap(h) )  
  return NArray.to_na([vor_dt*0, div_dt*0, h_dt])
}


# 出力ファイル名
ofn = ["vor.nc", "div.nc", "h.nc"]

gpna = NArray.to_na([vor.sp, div.sp, h.sp]) # GPhys変数をNArrayに格納


i = 0
while (dt*i < T) do
#while (i < 10000) do
  if (i%10 == 0) then # check CFL condition
    psi = invlap(gpna[0]) # 流線関数
    chi = invlap(gpna[1]) # 速度ポテンシャル
    u = (-psi.cos2_dy.gp + chi.dx.gp)/coslat # u
    v = (psi.dx.gp + chi.cos2_dy.gp)/coslat  # v
    umax = [u.abs.max.to_f, v.abs.max.to_f].max
    p "umax = #{umax}, CFL = #{(umax*dt/(dx*PI/180*R)).to_f}, #{i} steps. "
    auto_write(ofn[0], gpna[0].gp, i, false)
    auto_write(ofn[1], gpna[1].gp, i, false)
    auto_write(ofn[2], gpna[2].gp, i, false)
    
    binding.pry if File.exist?("stop")

  end

#    gpna = TimeIntegration.rk3(gpna,dt,1,&voreq2Dplane)
  gpna = TimeIntegration.rk3(gpna,dt,1,&shallow_sphere)

  p i*dt if (i%20 == 0)
  i += 1
  if gpna[0].real.max.to_f.nan? then
    raise "!!! Numerical Instability occurs !!!"
  end
    
end
