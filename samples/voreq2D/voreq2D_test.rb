require "numru/gphys"
require "numru/ganalysis"
require "../../lib/timeintegration.rb"
require "../../lib/gphysmath.rb"

require "pry"
include NumRu
include NMath
include GPMath

# 定数設定
dt = UNumeric[ARGV[0].to_f,"s"]   # 時間刻み幅
HD = UNumeric[5E4.to_f,"m2.s-1"]  # 拡散係数 (∇^2)
T = UNumeric[86400*5,"s"]         # 積分時間
Beta = UNumeric[0.0E-11, "m-1.s-1" ]
R = GAnalysis::Planet.radius
Omega = GAnalysis::Planet.omega

# 領域設定
# 周期境界の2次元平面
grid = GPhys.domain2D(
 {"range"=>0..1E7, "grid_points"=>64, "cyclic"=> true, "units"=>"m"},
 {"range"=>0..1E7, "grid_points"=>64, "cyclic"=> true, "units"=>"m"} )
# 2次元球面  
# grid = GPhys.domain2Dsphere(64,R) # 引数は経度方向の格子点数
  
# 変数設定（渦度）
vor = GPhys.make_var(grid,
  {"name"=>"vor", "long_name"=>"vorticity", "units"=>"s-1"} )

# 球面調和関数変換の準備
# GPhys.prepare_spht(vor)
# GPhys.use_spht

##############
# 初期分布設定 #
##############
# グリッドデータで与える場合
xc0 = 4000E3; yc0 = 5000E3; xc1 = 6000E3; yc1 = 5000E3; rlen = 1000E3
# xc0 = 150; yc0 = 0; xc1 = 210; yc1 = 0; rlen = 20
xval = vor.coord(0).val; yval = vor.coord(1).val
val = vor.val
grid.shape[1].times{|j|
  grid.shape[0].times{|i|
    r0 = sqrt( (xval[i]-xc0)**2 + (yval[j]-yc0)**2 )
    r1 = sqrt( (xval[i]-xc1)**2 + (yval[j]-yc1)**2 )
    val[i,j,0] = (exp((rlen-r0)/rlen*2) + exp((rlen-r1)/rlen*2))*1.0E-5
  }
}
vor.replace_val(val)

GPhys.prepare_fft_wavenumber(vor)

# Vorticity equation on double cyclic 2D-plane
voreq2Dplane = lambda {|gpna,t|
  GPhys.use_fft
  zeta = gpna[0]
  psi = invlap(zeta)
  zeta_tendency = -jacobian(psi,zeta) - Beta*psi.dx + HD*lap(zeta)
  return NArray.to_na([zeta_tendency])
}

# Vorticity equation on 2D-sphere
voreq2Dsphere = lambda {|gpna,t|
  zeta = gpna[0]
  psi = invlap(zeta)
  zeta_tendency = -jacobian(psi,zeta) - 2*Omega/R*psi.dx + HD*(lap(zeta) + 2/(R*R)*zeta)
  return NArray.to_na([zeta_tendency])
}

ofn1 = "vor.nc"
gpna = NArray.to_na([vor])

i = 0
# 時間発展ループ
while (dt*i < T) do
  auto_write(ofn1, gpna[0], i, false)
  gpna = TimeIntegration.rk3(gpna,dt,1,&voreq2Dplane)
#  gpna = TimeIntegration.rk3(gpna,dt,10,&voreq2Dsphere)
  p i*dt if (i%20 == 0)
  i += 1
  if gpna[0].max.to_f.nan? then
    raise "!!! Numerical Instability occurs !!!"
  end
end
