###
=begin
= NAME
TimeIntegration
= AUTHOR
Hiroki Kashimura
= DESCRIPTION
Time-Integration Methods for GPhys
= VERSION

=end

require "numru/gphys"

module TimeIntegration
  class << self

		# Arguments for the methods below
	  # gpna: NArray of GPhys objects of the data
	  # dt: time-step
	  # n (=1): number of steps to proceed
	  # &f: function defined in lambda

	  # euler (forward) scheme method
	  def euler(gpna,dt,n=1,&f)
	    t = gpna[0].coord(-1) # 時間軸(長さ1)を取り出す
	    var_num = gpna.length
	  	n.times{
		    k1na = f.call(gpna,t)*dt # 各変数の増分の計算
			  gpna = gpna + k1na # 次ステップの値
 	  		t  = t + dt
		  }
		  var_num.times{|v|  gpna[v].axis(-1).set_pos(t)  }
	    return gpna
	  end

	  # 4-stage Runge-Kutta
		def rk4(gpna,dt,n=1,&f)
	  	t = gpna[0].coord(-1) # 時間軸(長さ1)を取り出す
	  	var_num = gpna.length
	  	n.times{
			  k1na = f.call(gpna,t)*dt
			  k2na = f.call(gpna+k1na/2,t+dt/2)*dt
			  k3na = f.call(gpna+k2na/2,t+dt/2)*dt
			  k4na = f.call(gpna+k3na,t+dt)*dt
			  gpna = gpna + (k1na+k2na*2+k3na*2+k4na)/6
		  	t  = t + dt
			}
 		  var_num.times{|v|  gpna[v].axis(-1).set_pos(t)  }
		  return gpna
		end

		# 3-stage Runge-Kutta
		def rk3(gpna,dt,n=1,&f)
	  	t = gpna[0].coord(-1) # 時間軸(長さ1)を取り出す
	  	var_num = gpna.length
	  	n.times{
			  k1na = f.call(gpna,t)*dt
			  k2na = f.call(gpna+k1na/3,t+dt/3)*dt
			  k3na = f.call(gpna+k2na*2/3,t+dt*2/3)*dt
			  gpna = gpna + (k1na+k3na*3)/4
			  t  = t + dt
			}
			var_num.times{|v|  gpna[v].axis(-1).set_pos(t)  }
 		  return gpna
		end

		# 2-stage Runge-Kutta (Heun method)
		def rk2(gpna,dt,n=1,&f)
		  t = gpna[0].coord(-1) # 時間軸(長さ1)を取り出す
		  var_num = gpna.length
		  n.times{
			  k1na = f.call(gpna,t)*dt
			  k2na = f.call(gpna+k1na,t+dt)*dt
			  gpna = gpna + (k1na+k2na)/2
			  t  = t + dt
			}
			var_num.times{|v|  gpna[v].axis(-1).set_pos(t)  }
		  return gpna
		end

		# Leap-Frog Method
		def leapfrog(gpna0,gpna1,dt,n=1,filter=nil,&f)
		  t0 = gpna0[0].coord(-1) ; t1 = gpna1[0].coord(-1) # 時間軸(長さ1)を取り出す
		  var_num = gpna0.length
		  n.times{
		    gpna2 = gpna0 + f.call(gpna1,t1)*dt*2
		    if (filter == "asselin") then
		    	gpna1 = asselin(gpna0, gpna1, gpna2)
		    elsif (filter == "williams") then
		    	gpna1, gpna2 = williams(gpna0, gpna1, gpna2)
		    elsif (filter == nil ) then
		    else raise "value of 'filter' is invalid."
		    end
		    t0  = t0 + dt; t1  = t1 + dt
		    gpna0 = gpna1    ; gpna1 = gpna2
		  }
		  var_num.times{|v|
		    gpna0[v].axis(-1).set_pos(t0)
		    gpna1[v].axis(-1).set_pos(t1)
		  }
      return gpna0, gpna1
		end

		# trapezoidal method (Kurihara, 1965)
		def trapezoidal(gpna0,gpna1,dt,n=1,filter=nil,&f)
			t0 = gpna0[0].coord(-1) ; t1 = gpna1[0].coord(-1) # 時間軸(長さ1)を取り出す
		  var_num = gpna0.length
		  n.times{
		  	gpna_tmp = gpna0 + f.call(gpna1,t1)*dt*2
		    gpna2 = gpna0*0 + gpna1 + (f.call(gpna1, t1) + f.call(gpna_tmp, t1+dt))*dt/2 # なぜか最初の項がないと時刻がおかしくなる
		    if (filter == "asselin") then
		    	gpna1 = asselin(gpna0, gpna1, gpna2)
		    elsif (filter == "williams") then
		    	gpna1, gpna2 = williams(gpna0, gpna1, gpna2)
		    elsif (filter == nil ) then
		    else raise "value of 'filter' is invalid."
		    end
		    t0  = t0 + dt; t1  = t1 + dt
		    gpna0 = gpna1    ; gpna1 = gpna2
		  }
		  var_num.times{|v|
		    gpna0[v].axis(-1).set_pos(t0)
		    gpna1[v].axis(-1).set_pos(t1)
		  }
      return gpna0, gpna1
    end


		# Asselin filter
		def asselin(gpna0,gpna1,gpna2,delta=0.05)
		  gpna1 = (1.0-2.0*delta)*gpna1 +delta*(gpna2 + gpna0)
		  return gpna1
		end
		# Williams (2009) filter (RAW filter) # Williams (2009)論文とdeltaが２倍異なる
		def williams(gpna0,gpna1,gpna2,delta=0.05,alpha=0.5)
    	alpha = 0.5; delta = 0.05
		  dna = delta*(gpna0 - 2.0*gpna1 + gpna2)
	    gpna1 = gpna1 + dna*alpha
	    gpna2 = gpna2 - dna*(1.0-alpha)
			return gpna1, gpna2
		end

		# Adams-Bashforth 3rd-order
		def ab3(gpna0,gpna1,gpna2,dt,n=1,&f)
			t0 = gpna0[0].coord(-1) ; t1 = gpna1[0].coord(-1); t2 = gpna2[0].coord(-1)# 時間軸(長さ1)を取り出す
  	  var_num = gpna0.length
		  n.times{
			  gpna3 = gpna0*0 + gpna2 + (23.0*f.call(gpna2,t2) -16.0*f.call(gpna1,t1) + 5.0*f.call(gpna0,t0))*dt/12
			  t0 = t0 + dt ; t1 = t1 + dt ; t2 = t2 + dt
			  gpna0 = gpna1; gpna1 = gpna2; gpna2 = gpna3
			}
			var_num.times{|v|
		    gpna0[v].axis(-1).set_pos(t0)
		    gpna1[v].axis(-1).set_pos(t1)
		    gpna2[v].axis(-1).set_pos(t2)
		  }
		  return gpna0, gpna1, gpna2
		end

	end
end
