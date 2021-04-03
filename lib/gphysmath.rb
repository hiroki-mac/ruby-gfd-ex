###
=begin
= NAME
GPhysMath
= AUTHOR
Hiroki Kashimura
= DESCRIPTION
Mathematical methods for GPhys
= VERSION

=end

require "numru/gphys"
require "numru/ganalysis"
require "../../lib/spherical_harmonics_next.rb"
include NumRu

# 複数のGPhysデータを必要とする演算（ヤコビアンなど）はGPMathモジュールのメソッドとして定義する。
# 単一のGPhysデータに対する演算（d/dx, ラプラシアンなど）は、GPhysのクラスメソッドとして追加する。
# このうち、数学表記で変数の前に来るもの（ラプラシアンなど）は、GPMathモジュールのメソッドとしても定義する。
module GPMath
  @@use_fft = false
  @@use_spht = false
  @@output_in_spectral = false
  @@radius = GAnalysis::Planet.radius
  @@boundary_condition = [2,2,2] # 1:linear 2:cyclic
  @@x = 0; @@y = 1; @@z = 2
  @@nmax = nil
  @@fft_wavenumber_x = nil
  @@fft_wavenumber_y = nil
  @@fft_wavenumber_z = nil
  @@fft_wavenumber_for_lap = nil
  @@fft_wavenumber_for_invlap = nil
  @@inv_cos2lat = nil

  module_function
    def lap(gp)   ; gp.lap   ; end
    def invlap(gp); gp.invlap; end

    def jacobian(gp0,gp1)
      name0 = gp0.name; name1 = gp1.name
      gpjacob = gp0.dx*gp1.dy - gp0.dy*gp1.dx
      gpjacob.rename("J"+name0+name1)
      gpjacob.set_att("long_name", "J("+name0+","+name1+")")
      return gpjacob
    end

    # netcdfへの逐次書き出し
    def auto_write(outfilename, gary, index, compress=false) # 出力ファイル名, Array of GPhys, 時刻のindex, 圧縮の可否
      NetCDF.creation_format=(NetCDF::NC_NETCDF4 | NetCDF::NC_CLASSIC_MODEL) if compress

      gary = [gary] if gary.class == GPhys
      udname = gary[0].axnames[-1] # UNLIMITED次元の名前
      udval = gary[0].coordinate(udname).val

      if (index == 0) then # 初期設定
        outncfile = NetCDF.create(outfilename)
        gary[0].axnames.each{|n|
          a = gary[0].coordinate(n)
          if (n == udname) then
            outncfile.def_dim(n, 0) # 次元定義 [UNLIMITED次元]
          else
            outncfile.def_dim(n, a.length) # 次元定義
          end
          outncfile.def_var(n, a.ntype, [n]) # 次元変数定義
          a.att_names.each{|m| outncfile.var(n).put_att(m, a.get_att(m), nil) } # 次元変数属性の設定
        }
        gary.each{|g|
          outncfile.def_var(g.name, g.ntype, g.axnames) # 変数定義
          outncfile.var(g.name).deflate(1, true) if compress # set commpression level and shuffle.
          g.att_names.each{|n| outncfile.var(g.name).put_att(n, g.get_att(n), nil) } # 変数属性の設定
        }

        # set global atts if used in gpv.rb
        if (@sources) then
          if (@sources.length == 1 && @sources[0].include?(".nc")) then # copy the global atts of the source file.
            old_nc = NetCDF.open(@sources[0], "r")
            old_nc.each_att{|at| nc.put_att(at.name, at.get) }
          end
          outncfile.put_att("sources", @sources.join(", "))
          outncfile.put_att("command", @commandline )
          outncfile.put_att("working_directory", @current_dir )
          outncfile.put_att("process_date", "#{Time.now}" )
        end

        outncfile.enddef # 定義モード終了
        gary[0].axnames[0..-2].each{|n|
          outncfile.var(n).put(gary[0].coordinate(n).val) # 次元変数の値代入
        }
        outncfile.close # 一旦閉じる
      end

      outncfile = NetCDF.open(outfilename, "a") # 追記モードで開く
      st = outncfile.dim(udname).length
      et = outncfile.dim(udname).length + (udval.size - 1)
      outncfile.var(udname).put(udval,  "start"=>[st], "end"=>[et]) # UNLIMITED次元の値代入

      gary.each{|g|
        sta = Array.new(g.rank-1,0) << st
        eta = Array.new(g.rank-1,-1) << et
        outncfile.var(g.name).put(g.val, "start"=>sta,"end"=>eta) # 変数の値代入
      }
      outncfile.close # 処理の途中でも、別プロセスで可視化できるように毎回閉じる。
    end



end

################################################################################
################################################################################

class GPhys
  include GPMath
  include SphericalHarmonics

  def xdim; @@x; end
  def ydim; @@y; end
  def zdim; @@z; end

  def self.use_fft(bool=true); @@use_fft = bool; end
  def self.use_spht(bool=true); @@use_spht = bool; end
  def self.linear_bc(dim=nil)
    if dim == nil then @@boundary_condition = [1,1,1]
    else @@boundary_condition[dim] = 1 end
  end
  def self.cyclic_bc(dim=nil)
    if dim == nil then @@boundary_condition = [2,2,2]
    else @@boundary_condition[dim] = 2 end
  end
  def self.output_in_spectral(bool=true); @@output_in_spectral = bool; end

  def deriv(x) # 微分
    if @@use_spht then
      if (x == @@x) then # d/dlon
#        dgp = SphericalHarmonics.sh_dlon(self.gp)/@@radius
        smn = self.sp
        dgp = SphericalHarmonics.sp_dlon(smn)/@@radius # d/dlon * (1/radius)
      elsif (x == @@y) then # d/dmu * (1/radius) ; mu = sin(lat)
        if (self.ntype == "complex")
          dgp = SphericalHarmonics.sh_invtrans(self, false, true, false, self.size[1]-1)
        else
          smn = (self*@@inv_cos2lat).sp(@@nmax+1)
          dgp = SphericalHarmonics.sp_dmu_cos2(smn)/@@radius
        end
      end
      if @@output_in_spectral then
        return dgp.sp
      else
        return dgp.gp
      end
    elsif @@use_fft then  # TODO: 非周期境界のときでもFFTで計算可能にする
      self.fft_deriv(x)
    else
      self.threepoint_O2nd_deriv(x,@@boundary_condition[x])
    end
  end

  # 微分のエイリアス
  def dx; self.deriv(@@x); end
  def dy; self.deriv(@@y); end
  def dz; self.deriv(@@z); end

  def cos2_dy
    smn = self.sp(@@nmax+1)
    d_smn = SphericalHarmonics.sp_cos2_dmu(self)/@@radius
    if @@output_in_spectral then
      return d_smn
    else
      return d_smn.gp
    end
  end

  def get_coord_as_gpdata(d)
    gpdata = self.real.copy; val = gpdata.val
    coord = self.coord(d).val
    if (d == 0) then
      val = val*0 + coord
    elsif (d == 1) then
      val.shape[0].times{|i|
        val[i,false] = val[i,false]*0 + coord
      }
    elsif (d == 2) then
      val.shape[1].times{|j|
        val.shape[0].times{|i|
          val[i,j,false] = val[i,j,false]*0 + coord
        }
      }
    else
      raise "d > 2 is not supported."
    end
    gpdata.replace_val(val)
    gpdata.set_att("units",self.coord(d).units.to_s)
    gpdata.rename(self.coord(d).name)
    return gpdata
  end

  # get NArray of wavenumber from -K to K.
  def fft_wavenumber(d)
    sp = self.fft(false,d)
    kdata = sp.get_coord_as_gpdata(d)
    kval = kdata.val
    kmax = sp.coord(d).val[-1] + sp.coord(d).val[1]
    kval.map!{|k|
      if (k > kmax/2) then
        k - kmax
      else
        k
      end
    }
#    kdata.replace_val(kval)
    return kval
  end

  # prepare wavenumber NArray used in deriv, lap, invlap operations.
  def self.prepare_fft_wavenumber(gp)
    # 非周期境界なら鏡像拡張する
    gp = gp.mirror_ext(@@x) if (@@boundary_condition[@@x] == 1)
    gp = gp.mirror_ext(@@y) if (@@boundary_condition[@@y] == 1)
    gp = gp.mirror_ext(@@z) if (@@boundary_condition[@@z] == 1)

    rank = gp.rank
    @@fft_wavenumber_x = gp.fft_wavenumber(@@x) if (gp.shape[@@x] > 1)
    if (rank > @@y) then
      @@fft_wavenumber_y = gp.fft_wavenumber(@@y) if (gp.shape[@@y] > 1)
    end
    if (rank > @@z) then
      @@fft_wavenumber_z = gp.fft_wavenumber(@@z) if (gp.shape[@@z] > 1)
    end
    if (@@fft_wavenumber_x and @@fft_wavenumber_y) then
      @@fft_wavenumber_for_lap = -(@@fft_wavenumber_x*@@fft_wavenumber_x + @@fft_wavenumber_y*@@fft_wavenumber_y) # -(k*k + l*l)
      @@fft_wavenumber_for_invlap = -1.0/(@@fft_wavenumber_for_lap + @@fft_wavenumber_for_lap.eq(0)) # -1/(k*k + l*l)
    end
  end

  def self.prepare_spht(gp)
    lon = gp.coord(@@x).val
    lat = gp.coord(@@y).val
    @@nmax = (lon.length-1)/3
    SphericalHarmonics.sh_init(lon,lat,@@nmax,gp)
    coslat = SphericalHarmonics.coslat()
    @@inv_cos2lat = 1.0/(coslat*coslat)
    return true
  end

  # derivative by FFT
  # d: number of the dimension for derivative
  # n: 微分の階数
  def fft_deriv(d=0,n=1)
    gp = self.copy
    gp = gp.mirror_ext(@@x) if (@@boundary_condition[@@x] == 1)
    gp = gp.mirror_ext(@@y) if (@@boundary_condition[@@y] == 1)
    gp = gp.mirror_ext(@@z) if (@@boundary_condition[@@z] == 1)
    sp = gp.fft(false,d)
    case d
    when @@x;      kval = @@fft_wavenumber_x
    when @@y;      kval = @@fft_wavenumber_y
    when @@z;      kval = @@fft_wavenumber_z
    else;          raise "not supported."
    end
    coef = kval.to_type("complex")*Complex::I
    sp = sp*coef**n
    gpout = sp.fft(true,d).real
    # <<set attribute>>
    axis = self.coord(d)
    gpout.rename("d#{self.name}_d#{axis.name}")       # ex. "dT_dx"
    gpout.units = self.units/axis.units        # set units
    if self.get_att('long_name') && axis.get_att('long_name') then
      long_name = "d_#{self.get_att('long_name')}_d_#{axis.get_att('long_name')}"
    else
      long_name = name
    end
    gpout.set_att("long_name",long_name)            # set long_name
    
    gpout = gpout.cut_mirror_ext(@@x, self) if (@@boundary_condition[@@x] == 1)
    gpout = gpout.cut_mirror_ext(@@y, self) if (@@boundary_condition[@@y] == 1)
    gpout = gpout.cut_mirror_ext(@@z, self) if (@@boundary_condition[@@z] == 1)

    return gpout
  end

  # Laplacian on x-y plane
  def lap
    name = self.name
    if (@@use_spht) then
      smn = self.sp
      gplap = SphericalHarmonics.sp_laplacian(smn)/(@@radius*@@radius) # スペクトル→スペクトル
      gplap = gplap.gp unless @@output_in_spectral

    elsif (@@use_fft) then
      if (@@boundary_condition[@@x] == 2 && @@boundary_condition[@@y] == 2) then # 2重周期
        m_ext = self
      elsif (@@boundary_condition[@@x] == 1 && @@boundary_condition[@@y] == 1) then
        m_ext = self.mirror_ext(@@x).mirror_ext(@@y) # x, y ともに非周期
      elsif (@@boundary_condition[@@x] == 1) then
        m_ext = self.mirror_ext(@@x) # xは非周期 y は周期
      else
        m_ext = self.mirror_ext(@@y) # xは周期 y は非周期
      end
      sp = m_ext.fft(false,@@x,@@y)
      sp = sp.mul!(@@fft_wavenumber_for_lap)
      gplap = sp.fft(true,@@x,@@y).real
      gplap.units = self.units/self.coord(0).units**2
      gplap = gplap.cut_mirror_ext(@@x, self) if (@@boundary_condition[@@x] == 1)
      gplap = gplap.cut_mirror_ext(@@y, self) if (@@boundary_condition[@@y] == 1)
    else
      gplap = self.deriv2nd(@@x,@@boundary_condition[@@x]) + self.deriv2nd(@@y,@@boundary_condition[@@y])
    end

    gplap.rename("lap_"+name)
    gplap.set_att("long_name","lap_"+name)
    return gplap
  end

  # inverse of Laplacian on x-y plane
  def invlap
    name = self.name
    if (@@use_spht) then # 球面
      if (self.ntype == "complex") then
        gpout = SphericalHarmonics.sp_invlaplacian(self)*(@@radius*@@radius) # スペクトル→スペクトル
      else
        gpout = SphericalHarmonics.sh_invlaplacian(self)*(@@radius*@@radius) # グリッド→スペクトル→グリッド
      end
    elsif (@@use_fft) then
      if (@@boundary_condition[@@x] == 2 && @@boundary_condition[@@y] == 2) then # 2重周期
        m_ext = self
      elsif (@@boundary_condition[@@x] == 1 && @@boundary_condition[@@y] == 1) then
        m_ext = self.mirror_ext(@@x).mirror_ext(@@y) # x, y ともに非周期
      elsif (@@boundary_condition[@@x] == 1) then
        m_ext = self.mirror_ext(@@x) # xは非周期 y は周期
      else
        m_ext = self.mirror_ext(@@y) # xは周期 y は非周期
      end
      sp = m_ext.fft(false,@@x,@@y)
      sp = sp.mul!(@@fft_wavenumber_for_invlap)
      gpout = -sp.fft(true,@@x,@@y).real
      gpout.units = self.units*self.coord(0).units**2
      gpout = gpout.cut_mirror_ext(@@x, self) if (@@boundary_condition[@@x] == 1)
      gpout = gpout.cut_mirror_ext(@@y, self) if (@@boundary_condition[@@y] == 1)
    else # 差分法
      # gpout = self.int2(@@x).int2(@@y) これは間違い
      # SOR法 境界値が固定になるので注意（周期境界には非対応）
      require "numru/ganalysis/pde"
      gpout = self.copy
      f = self.val.reshape(self.shape[0], self.shape[1])
      z = f.clone
      a = z*0.0 + 1.0
      b = z*0.0
      dx = self.coord(@@x).val[1] - self.coord(@@x).val[0]
      dy = self.coord(@@y).val[1] - self.coord(@@y).val[0]
      res = 0.0
      GAnalysis::PDE.SOR_2D_2ndorder(z,a,b,a,b,b,f,dx,dy,1.9,eps:1E-8)
      gpout.replace_val(z.reshape(*(self.shape)))
      gpout.units = self.units*self.coord(@@x).units*self.coord(@@y).units
    end
    gpout.rename("invlap_"+name)
    gpout.set_att("long_name","invlap_"+name)
  end

  # spherical harmics spectral transform (grid => spectral)
  def sp(nmax=nil)
    case self.ntype
    when "complex"
      return self
    else
      if nmax == nil then
        sp = SphericalHarmonics.sh_trans(self)
      else
        sp = SphericalHarmonics.sh_trans(self,false,false,nmax)
      end
      return sp
    end
  end

  # spherical harmonics spectral transform (spectral => grid)
  def gp(nmax=nil)
    case self.ntype
    when "complex"
      nmax = self.shape[1] - 1 if nmax == nil
      gp = SphericalHarmonics.sh_invtrans(self,false,false,false,nmax)
      gp.rename(gp.name.delete("_sp_g"))
      return gp
    else
      return self
    end
  end

  # Tools for making Axis and Grid objects easily.
  def self.makeAxis(opts={},na=nil)
    default_opts = {
      "range" => 0..1,
      "grid_points" => 64,
      "units" => "1",
      "name" => "x",
      "long_name" => "x-axis",
      "cyclic" => true,
      "virtual_points" => 0,
      "cell_center" => true
    }
    default_opts.each_key{|key|
      opts[key] = default_opts[key] if opts[key].nil?
    }
    st = opts["range"].begin; ed = opts["range"].end
    n = opts["grid_points"]
    vn = opts["virtual_points"]
    if (opts["cyclic"] or n == 1) then
      dx = (ed - st).to_f/n
    else
      dx = (ed - st).to_f/(n-1)
    end
    if (na) then
      x_na = na
    else
      if (opts["cyclic"]) then
        x_na = NArray.float(n+vn*2).indgen!*dx + st - vn*dx
      elsif (opts["cell_center"]) then
        x_na = NArray.float(n+vn*2).indgen!*dx + st - vn*dx + dx*0.5
      else # cell boundary
        x_na = NArray.float(n+vn*2).indgen!*dx + st - vn*dx
      end
    end

    opts.each_key{|key|
      unless ([String, Array, NArray].include?(opts[key].class)) then
        opts[key] = opts[key].to_s
      end
    }
    x_va = VArray.new(x_na, opts, opts["name"])
    return Axis.new.set_pos(x_va)
  end

  def self.domain1D(opts={})
    @@boundary_condition[@@x] = 1 if (xopts["cyclic"] == false)
    xaxis = GPhys.makeAxis(opts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})
    return Grid.new(xaxis,taxis)
  end

  def self.domain2D(xopts={},yopts={})
    @@boundary_condition[@@x] = 1 if (xopts["cyclic"] == false)
    @@boundary_condition[@@y] = 1 if (yopts["cyclic"] == false)
    xaxis = GPhys.makeAxis(xopts)
    default_yopts = {"name" => "y", "long_name" => "y-axis"}
    default_yopts.each_key{|key|
      yopts[key] = default_yopts[key] if yopts[key].nil?
    }
    yaxis = GPhys.makeAxis(yopts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})
    return Grid.new(xaxis,yaxis,taxis)
  end

  def self.domain3D(xopts={},yopts={},zopts={})
    @@boundary_condition[@@x] = 1 if (xopts["cyclic"] == false)
    @@boundary_condition[@@y] = 1 if (yopts["cyclic"] == false)
    @@boundary_condition[@@z] = 1 if (zopts["cyclic"] == false)
    xaxis = GPhys.makeAxis(xopts)
    default_yopts = {"name" => "y", "long_name" => "y-axis"}
    default_yopts.each_key{|key|
      yopts[key] = default_yopts[key] if yopts[key].nil?
    }
    yaxis = GPhys.makeAxis(yopts)
    default_zopts = {"name" => "z", "long_name" => "z-axis"}
    default_zopts.each_key{|key|
      zopts[key] = default_zopts[key] if zopts[key].nil?
    }
    zaxis = GPhys.makeAxis(zopts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})
    return Grid.new(xaxis,yaxis,zaxis,taxis)
  end

  def self.domain2Dsphere(imax,radius=nil)
    xaxis = GPhys.makeAxis({"range"=>0..360, "grid_points"=>imax, "name"=>"lon", "units"=>"deg", "long_name"=>"longitude", "cyclic"=>true})
    mu, weight = SphericalHarmonics.gauss_w(imax/2)
    gauss_lat = SphericalHarmonics.mu2deg(mu)
    yaxis = GPhys.makeAxis({"range"=>-90..90, "grid_points"=>imax/2, "name"=>"lat", "units"=>"deg", "long_name"=>"latitude", "cyclic"=>false}, gauss_lat)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})

    if (radius) then
      @@radius = radius
    end

    return Grid.new(xaxis,yaxis,taxis)
  end

  def self.domain3Dsphere(imax,zopts={},radius=nil)
    xaxis = GPhys.makeAxis({"range"=>0..360, "grid_points"=>imax, "name"=>"lon", "units"=>"deg", "long_name"=>"longitude", "cyclic"=>true})
    mu, weight = SphericalHarmonics.gauss_w(imax/2)
    gauss_lat = SphericalHarmonics.mu2deg(mu)
    yaxis = GPhys.makeAxis({"range"=>-90..90, "grid_points"=>imax/2, "name"=>"lat", "units"=>"deg", "long_name"=>"latitude", "cyclic"=>false}, gauss_lat)

    @@boundary_condition[@@z] = 1
    default_zopts = {"name" => "z", "long_name" => "height"}
    default_zopts.each_key{|key|
      zopts[key] = default_zopts[key] if zopts[key].nil?
    }
    zaxis = GPhys.makeAxis(zopts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})

    if (radius) then
      @@radius = radius
    end

    return Grid.new(xaxis,yaxis,zaxis,taxis)
  end


  def self.make_var(grid,opts={})
    default_opts = {
      "units" => "1",
      "name" => "var",
      "long_name" => "var_name"
    }
    default_opts.each_key{|key|
      opts[key] = default_opts[key] if opts[key].nil?
    }
    na = NArray.float(*grid.shape)
    va = VArray.new(na, opts, opts["name"])
    return GPhys.new(grid,va)
  end

  def mirror_ext(dim) # GPhys オブジェクトの 鏡像拡張 (正負反転)
    mirror = self.copy
    mirror.axis(dim).set_pos(-mirror.coord(dim)) # 軸の正負反転
    mirror.coord(dim).replace_val(mirror.coord(dim).val.reverse) # 軸の値(順序)を反転
    mirror.replace_val(-mirror.val.reverse(dim)) # 変数の値を軸に当該沿って反転し、正負も反転させる
    if (mirror.coord(dim).val[-1] == self.coord(dim).val[0]) # 0の重複を除く
      tmp = mirror.cut_rank_conserving({self.axis(dim).name=>mirror.coord(dim).val[0]})
      dx = self.coord(dim)[1].val - self.coord(dim)[0].val
      tmp.coord(dim).replace_val(tmp.coord(dim).val - dx)
      tmp.replace_val(tmp.val*0.0)
      return GPhys.join([tmp, mirror[true,0..-2,false], self])
    else
      return GPhys.join([mirror, self])
    end
  end
  
  def cut_mirror_ext(dim, gp_org) # 鏡像拡張した部分をカットする。gp_org は元々のGPhysオブジェクト
    return self.cut({gp_org.axis(dim).name => self.coord(dim).val[-gp_org.shape[dim]]..self.coord(dim).val[-1]})
  end

  def int(dim, bnd_v=nil)

    gp = self.copy
    axis = self.coord(dim).val; im = axis.length
    flag_mean = false
    if (bnd_v == nil) then; bnd_v = 0.0; flag_mean = true; end
    
    case dim
    when 0
      gp[0,false] = gp[0,false]*0.0 + bnd_v # 下端の値
      gp[1,false] = gp[0,false] + (axis[1]-axis[0])*self[0,false]
      for i in 1..(im-2)
        gp[i+1,false] = gp[i-1,false] + (axis[i+1]-axis[i-1])*self[i,false]
      end
    when 1
      gp[true,0,false] = gp[true,0,false]*0.0 + bnd_v # 下端の値
      gp[true,1,false] = gp[true,0,false] + (axis[1]-axis[0])*self[true,0,false]
      for i in 1..(im-2)
        gp[true,i+1,false] = gp[true,i-1,false] + (axis[i+1]-axis[i-1])*self[true,i,false]
      end
    when 2
      gp[true,true,0,false] = gp[true,true,0,false]*0.0 + bnd_v # 下端の値
      gp[true,true,1,false] = gp[true,true,0,false] + (axis[1]-axis[0])*self[true,true,0,false]
      for i in 1..(im-2)
        gp[true,true,i+1,false] = gp[true,true,i-1,false] + (axis[i+1]-axis[i-1])*self[true,true,i,false]
      end
    else; raise "not supported."
    end
    gp.units=(gp.units * gp.coord(dim).units)
    gp = gp.eddy(dim) if flag_mean
    return gp
  
  end

  def intO2(dim, bnd_v=nil)
    gp = self.copy
    axis = self.coord(dim).val; im = axis.length
    flag_mean = false
    if (bnd_v == nil) then; bnd_v = 0.0; flag_mean = true; end
    case dim
    when 0
      gp[0,false] = gp[0,false]*0.0 + bnd_v # 下端の値
      gp[1,false] = (axis[1]-axis[0])*self[0,false]
      for i in 1..(im-2)
        # t = x_{i+1} - x_{i}, s = x_{i} - x_{i-1}
        tbs = (axis[i+1]-axis[i])/(axis[i]-axis[i-1]) # t/s
        spt = axis[i+1] - axis[i-1] # s + t
        gp[i+1,false] = (1.0-tbs*tbs)*gp[i,false] + (tbs*tbs)*gp[i-1,false] + (tbs*spt)*self[i,false]
      end
    when 1
    when 2
    else; raise "not supported."
    end
    gp.units=(gp.units * gp.coord(dim).units)
    gp = gp.eddy(dim) if flag_mean
    return gp
  end

  # 2階積分
  def int2(dim, sb=nil, eb=nil) # sbnd 始点境界の値、 ebnd 終点境界の値
    gp = self.copy
    axis = self.coord(dim).val; im = axis.length; flag_mean = false
    if (sb == nil) then; sb = 0.0; flag_mean = true; end
    case @@boundary_condition[dim]
    when 1 # 非周期境界条件
      eb = sb if (eb == nil)
    when 2 # 周期境界条件
      eb = sb
    end
    
    case dim
    when 0
      gp[0,false] = gp[0,false]*0.0 + sb # 下端の値
      gp[1,false] = gp[0,false] + 0.5*(axis[1]-axis[0])*(axis[1]-axis[0])*self[0,false] # gp_(-1) = gp_1 を仮定
      for i in 1..(im-2)
        # t = x_{i+1} - x_{i}, s = x_{i} - x_{i-1}
        tbs = (axis[i+1]-axis[i])/(axis[i]-axis[i-1]) # t/s
        sptt = (axis[i+1] - axis[i-1])*(axis[i+1]-axis[i])*0.5 # (s + t)*t/2
        gp[i+1,false] = (1.0+tbs)*gp[i,false] - (tbs)*gp[i-1,false] + (sptt)*self[i,false]
      end
      for i in 0..(im-1) # 積分定数の C*x 分の調整
        gp[i,false] = gp[i,false] - (gp[im-1,false] - gp[0,false] - eb)*(axis[i]-axis[0])/(axis[im-1]-axis[0])
      end
    when 1
      gp[true,0,false] = gp[true,0,false]*0.0 + sb # 下端の値
      gp[true,1,false] = gp[true,0,false] + 0.5*(axis[1]-axis[0])*(axis[1]-axis[0])*self[true,0,false] # gp_(-1) = gp_1 を仮定
      for i in 1..(im-2)
        # t = x_{i+1} - x_{i}, s = x_{i} - x_{i-1}
        tbs = (axis[i+1]-axis[i])/(axis[i]-axis[i-1]) # t/s
        sptt = (axis[i+1] - axis[i-1])*(axis[i+1]-axis[i])*0.5 # (s + t)*t/2
        gp[true,i+1,false] = (1.0+tbs)*gp[true,i,false] - (tbs)*gp[true,i-1,false] + (sptt)*self[true,i,false]
      end
      for i in 0..(im-1) # 積分定数の C*x 分の調整
        gp[true,i,false] = gp[true,i,false] - (gp[true,im-1,false] - gp[true,0,false] - eb)*(axis[i]-axis[0])/(axis[im-1]-axis[0])
      end
    when 2
      gp[true,true,0,false] = gp[true,true,0,false]*0.0 + sb # 下端の値
      gp[true,true,1,false] = gp[true,true,0,false] + 0.5*(axis[1]-axis[0])*(axis[1]-axis[0])*self[true,true,0,false] # gp_(-1) = gp_1 を仮定
      for i in 1..(im-2)
        # t = x_{i+1} - x_{i}, s = x_{i} - x_{i-1}
        tbs = (axis[i+1]-axis[i])/(axis[i]-axis[i-1]) # t/s
        sptt = (axis[i+1] - axis[i-1])*(axis[i+1]-axis[i])*0.5 # (s + t)*t/2
        gp[true,true,i+1,false] = (1.0+tbs)*gp[true,true,i,false] - (tbs)*gp[true,true,i-1,false] + (sptt)*self[true,true,i,false]
      end
      for i in 0..(im-1) # 積分定数の C*x 分の調整
        gp[true,true,i,false] = gp[true,true,i,false] - (gp[true,true,im-1,false] - gp[true,true,0,false] - eb)*(axis[i]-axis[0])/(axis[im-1]-axis[0])
      end
    else; raise "not supported."
    end
    gp.units=(gp.units * gp.coord(dim).units * gp.coord(dim).units )
    gp = gp.eddy(dim) if flag_mean
    return gp
  end



end
