# 1次元移流方程式
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
C = UNumeric[100.0, "m.s-1"]
# 領域設定
# 1次元周期境界
grid = GPhys.domain1D(
 {"range"=>0..4E7, "grid_points"=>256, "cyclic"=> true, "units"=>"m"})
 
# 変数設定（渦度）
zeta = GPhys.make_var(grid,
  {"name"=>"zeta", "long_name"=>"zeta", "units"=>"1"} )

##############
# 初期分布設定 #
##############
# グリッドデータで与える場合
xc0 = 4000E3; xc1 = 6000E3; rlen = 1000E3
# xc0 = 150; yc0 = 0; xc1 = 210; yc1 = 0; rlen = 20
xval = zeta.coord(0).val; val = zeta.val
grid.shape[0].times{|i|
  r0 = sqrt( (xval[i]-xc0)**2 )
  val[i,0] = (exp((rlen-r0)/rlen*2) )*1.0E-5
}
zeta.replace_val(val)

GPhys.prepare_fft_wavenumber(zeta)

# advection equation on double cyclic 1D 
advect1D = lambda {|gpna,t|
  GPhys.use_fft
  zeta = gpna[0]
  zeta_tendency = - C * zeta.dx + HD*zeta.dx.dx
  return NArray.to_na([zeta_tendency])
}

ofn1 = "zeta.nc"
gpna = NArray.to_na([zeta])

i = 0
# 時間発展ループ
while (dt*i < T) do
  auto_write(ofn1, gpna[0], i, false)
  gpna = TimeIntegration.rk4(gpna,dt,1,&advect1D)
  p i*dt if (i%20 == 0)
  i += 1
  if gpna[0].max.to_f.nan? then
    raise "!!! Numerical Instability occurs !!!"
  end
end
