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
  @@boundary_condition = 2 # 1:linear 2:cyclic
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
  def self.linear_bc; @@boundary_condition = 1; end
  def self.cyclic_bc; @@boundary_condition = 2; end
  def self.output_in_spectral(bool=true); @@output_in_spectral = bool; end

  def deriv(x) # 微分
    if @@use_spht then
      if (x == @@x) then # d/dlon
#        binding.pry
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
    elsif @@use_fft then
      self.fft_deriv(x)
    else
      self.threepoint_O2nd_deriv(x,@@boundary_condition)
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
    sp = self.fft(false,d)
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
      sp = self.fft(false,@@x,@@y)
      sp = sp.mul!(@@fft_wavenumber_for_lap)
      gplap = sp.fft(true,@@x,@@y).real
      gplap.units = self.units/self.coord(0).units**2
    else
      gplap = self.deriv2nd(@@x,@@boundary_condition) + self.deriv2nd(@@y,@@boundary_condition)
    end

    gplap.rename("lap_"+name)
    gplap.set_att("long_name","lap_"+name)
    return gplap
  end

  # inverse of Laplacian on x-y plane
  def invlap
    name = self.name
    if (@@use_spht) then
      if (self.ntype == "complex") then
        gpout = SphericalHarmonics.sp_invlaplacian(self)*(@@radius*@@radius) # スペクトル→スペクトル
      else
        gpout = SphericalHarmonics.sh_invlaplacian(self)*(@@radius*@@radius) # グリッド→スペクトル→グリッド
      end
    else
      sp = self.fft(false,@@x,@@y)
      sp = sp.mul!(@@fft_wavenumber_for_invlap)
      gpout = sp.fft(true,@@x,@@y).real
      gpout.units = self.units*self.coord(0).units**2
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
      "cyclic" => true
    }
    default_opts.each_key{|key|
      opts[key] = default_opts[key] unless opts[key]
    }
    st = opts["range"].begin; ed = opts["range"].end
    n = opts["grid_points"]
    if (opts["cyclic"] or n == 1) then
      dx = (ed - st).to_f/n
    else
      dx = (ed - st).to_f/(n-1)
    end
    opts.each_key{|key|
      unless ([String, Array, NArray].include?(opts[key].class)) then
        opts[key] = opts[key].to_s
      end
    }
    if (na) then
      x_na = na
    else
      x_na = NArray.float(n).indgen!*dx + st
    end
    x_va = VArray.new(x_na, opts, opts["name"])
    return Axis.new.set_pos(x_va)
  end

  def self.domain1D(opts={})
    xaxis = GPhys.makeAxis(opts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})
    return Grid.new(xaxis,taxis)
  end

  def self.domain2D(xopts={},yopts={})
    xaxis = GPhys.makeAxis(xopts)
    default_yopts = {"name" => "y", "long_name" => "y-axis"}
    default_yopts.each_key{|key|
      yopts[key] = default_yopts[key] unless yopts[key]
    }
    yaxis = GPhys.makeAxis(yopts)
    taxis = GPhys.makeAxis({"range"=>0..0, "grid_points"=> 1, "name"=>"t", "units"=>"s", "long_name"=>"time", "cyclic"=>false})
    return Grid.new(xaxis,yaxis,taxis)
  end

  def self.domain3D(xopts={},yopts={},zopts={})
    xaxis = GPhys.makeAxis(xopts)
    default_yopts = {"name" => "y", "long_name" => "y-axis"}
    default_yopts.each_key{|key|
      yopts[key] = default_yopts[key] unless yopts[key]
    }
    yaxis = GPhys.makeAxis(yopts)
    default_zopts = {"name" => "z", "long_name" => "z-axis"}
    default_zopts.each_key{|key|
      zopts[key] = default_zopts[key] unless zopts[key]
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

    default_zopts = {"name" => "z", "long_name" => "height"}
    default_zopts.each_key{|key|
      zopts[key] = default_zopts[key] unless zopts[key]
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
      opts[key] = default_opts[key] unless opts[key]
    }
    na = NArray.float(*grid.shape)
    va = VArray.new(na, opts, opts["name"])
    return GPhys.new(grid,va)
  end








end
