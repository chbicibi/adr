#! ruby -Ku
require 'fileutils'

# COL = {0=>"red", 1=>"blue", 2=>"green", 3=>"violet", 4=>"purple", 5=>"cyan"}
COL = ["red", "blue", "green", "purple", "cyan", "violet"]

def setterm data, term
  if String === term
    data.sub! /^term.*$/, "term = 1"
    data.sub /^file_image.*$/, "file_image = '#{term}'"
  elsif term.to_i > 0
    data.sub! /^term.*$/, "term = 1"
  else
    data.sub /^term.*$/, "term = #{term}"
  end
end

def setifo data, ifile
  data.sub /^file_orbit.*$/, "file_orbit = '#{ifile}'"
end

def setifp data, ifile
  data.sub /^file_point.*$/, "file_point = '#{ifile}'"
end

def setscale data, scale
  data.sub /^scale.*$/, "scale = #{scale}"
end

def setview data, angle1, angle2
  angle = {"1" => angle1, "2" => angle2}
  data.gsub(/^angle([12]).*$/){"angle#{$1} = #{angle[$1]}"}
end

def setorbit data, a
  if a.size > COL.size
    od = a.map{|i| "replot file_orbit every :::#{i + 1}::#{i + 1} w l lw 3 lc rgb 'green'\n"}
  else
    od = a.map{|i| "replot file_orbit every :::#{i + 1}::#{i + 1} w l lw 5 lc rgb '#{COL[i]}'\n"}# % ((i >= 1 and i <= 2) ? " dt '.'" : (i == 3 ? " dt '-'" : ""))}
  end
  data.sub(/^(##+)$/){od.join + $1}
end

def setpoint data, a, vx, vz
  vp = [Math.sin(vx) * Math.sin(vz), -Math.sin(vx) * Math.cos(vz), Math.cos(vx)].map{|a| a * 100}
  # vp = [0, -Math.sin(vx), Math.sin(vx)].map{|a| a * 10000}

  if a.size > COL.size
    od = a.map{|i| (s = "file_point%s every ::#{i}::#{i}:0 w p pt 7 ps %.1f lc rgb '%s'") and "replot %s, %s\n" % [s % ["", 3, "black"], s % [" using ($1%+d):($2%+d):($3%+d)" % vp, 2, "green"]]}
  else
    od = a.map{|i| (s = "file_point using ($1%+d):($2%+d):($3%+d) every ::#{i}::#{i}:0 w p pt 7 ps %.1f lc rgb '%s'") and "replot %s, %s\n" % [s % [*vp, 5, "black"], s % [*vp.map{|a| a * 2}, 3, COL[i]]]}
  end
  data.sub(/^(##+)$/){od.join + $1}
end

def setarrow data, n, pos, vec
  as = "set style arrow #{n + 1} size screen 0.02, 20 filled lw 5 lc rgb '#{COL[n]}' front"
  od = "set arrow #{n + 1} from %s to %s as #{n + 1}" % [pos.map(&:to_s), pos.zip(vec).map{|a| a.inject(:+).to_s}].map{|a| a.join(", ")}
  data.sub(/^(##+)$/){"#{as}\n#{od}\n#{$1}"}
end

def setmanual data, txt
  data.sub(/^(##+)$/){"#{txt}\n#{$1}"}
end

def plot data, index = nil
  plt = "plot#{index}.plt"
  File.write plt, data
  system "wgnuplot #{plt}"
  # system "echo #{data} | gnuplot -p"
  File.delete plt if index
end

# ------------------------------------------------------------------------------

input      = $*.reject{|e| e =~ /^\//}
opts       = $* - input
time_start = Time.now

case input[0]
when "0"
  plot File.read "plot.plt"

when "1"
  plot File.read "plot_orbit.plt"

when "2"
  plot File.read "plot_dv.plt"


when "otest"
  idata = File.read "plot_orbit.plt"
  10.times do |i|
    angle = File.readlines("C:/Yamada/Studies/MyEC/test_%02d.out" % i)[0][1..-1].to_f
    data = setterm  idata, 1
    data = setifo   data, "C:/Yamada/Studies/MyEC/test_%02d.out" % i
    data = setview  data, 60, angle - 90
    data = setorbit data, [*0...5]
    # data = setarrow data, 1, [7078, 0, 0], [2000, 0, 0]
    plot data
    File.rename "plot.png", "orbit_%03d.png" % i
  end
  system "convert -resize 800 -layers optimize -loop 0 -delay 100 -dispose 2 -transparent white orbit_*.png orbit_anim.gif"

when "o"
  idata = File.read "plot_orbit.plt"
  10.times do |i|
    data = setterm  idata, 1
    data = setifo   data, "plot/orbit#{i}.dat"
    data = setview  data, 60, 45
    data = setorbit data, [*0..2]
    # data = setarrow data, 1, [7078, 0, 0], [2000, 0, 0]
    plot data
    File.rename "plot.png", "plot/orbit_%03d.png" % i
  end

when "o1"
  p (File.readlines "arrow.dat")[0].strip.split(/\s+/).map(&:to_f).each_slice(3).to_a
  # exit
  idata = File.read "plot_orbit.plt"
  data  = setterm  idata, (input[1] or 0)
  data  = setifo   data, "orbit.dat"
  data  = setview  data, 60, 120
  data  = setorbit data, [0, 1]
  data  = setpoint data, [1]
  (0..3).each{|i| data  = setarrow data, i, *(File.readlines "arrow.dat")[i].strip.split(/\s+/).map(&:to_f).each_slice(3)}
  plot data

when "o2"
  p (File.readlines "arrow.dat")[0].strip.split(/\s+/).map(&:to_f).each_slice(3).to_a
  # exit
  idata = File.read "plot_orbit.plt"
  data  = setterm  idata, (input[1] or 0)
  data  = setifo   data, "orbit.dat"
  data  = setview  data, 60, 180
  data  = setorbit data, [1]
  # data  = setpoint data, [1]
  [2, 4, 5].each{|i| data  = setarrow data, i, *(File.readlines "arrow.dat")[i].strip.split(/\s+/).map(&:to_f).each_slice(3)}
  plot data

when "o3"
  exit unless input[1]
  index = input[1].to_i
  threads = []
  queue   = Queue.new
  10.times{|i| queue.push 0}

  [*index..index].each do |n|
    idata = File.read "plot_orbit.plt"
    FileUtils.mkdir_p "Plot/o#{n}"
    1.times do |i|
      sleep 0.2
      threads << Thread.new(n, i) do |n, i|
        queue.pop
        p [n, i] if i % 10 == 0
        of   = "plot/o#{n}/orbit_%03d.png" % i
        data = setterm  idata, of
        data = setifo   data,  "plot/data#{n}/orbit#{i}.dat"
        data = setifp   data,  "plot/data#{n}/point#{i}.dat"
        # data = setscale data,  5.0    if n / 3 == 0
        # data = setscale data,  3.8    if n / 3 == 1
        # data = setscale data,  7.0    if n / 3 == 2
        data = setscale data,  7.0
        # data = setview  data,  0, 0   if n % 3 == 0
        # data = setview  data,  90, 0  if n % 3 == 1
        # data = setview  data,  60, 45 if n % 3 == 2
        data = setview  data,  60, 45
        data = setview  data,  0, 0
        data = setorbit data,  [*0..4]
        # data = setpoint data,  [*0..4], 0,            0            if n % 3 == 0
        # data = setpoint data,  [*0..4], Math::PI / 2, 0            if n % 3 == 1
        # data = setpoint data,  [*0..4], Math::PI / 3, Math::PI / 4 if n % 3 == 2
        # data = setmanual data, "set label 3 left at screen 0.85, 0.1 '{/Arial:Bold=20 day #{i}}' textcolor rgb 'black' front"
        # data.sub! /set term png.+$/, "set term postscript pdf enhanced size width, height"
        data = data.sub(/(^\s*set term png).+$/){"#{$1}cairo enhanced truecolor fontscale 3.0 size width, height"}
        plot data, n * 100 + i
        # system "convert%s -density 600 -resize 800 #{of} plot/#{n}/cnv_pdf_#{i}.png" % (i == 0 ? " -alpha off" : "")
        system "convert -resize 800 #{of} #{of}"
        queue.push 0
      end
    end
    threads.each{|t| t.join}
    # system "convert -density 1500 -geometry 800x600 -transparent white -dispose previous -layers optimize -loop 0 -delay 10 plot/#{n}/orbit_*.pdf plot/orbit_#{n}.gif"
    # system "convert -density 1500 -geometry 800x600 -dispose previous -layers optimize -loop 0 -delay 10 plot/#{n}/*.pdf plot/orbit_pdf_#{n}.gif"
    # system "convert -layers optimize -loop 0 -delay 10 -transparent white -dispose previous plot/#{n}/*.png plot/orbit_png_#{n}.gif"
    # system "convert -layers optimize -loop 0 -delay 10 plot/#{n}/*.png plot/orbit_png_#{n}.gif"
  end

when "a"
  exit unless input[1]
  threads = []
  queue   = Queue.new
  3.times{|i| queue.push 0}

  index = input[1].to_i

  [*index..index].each do |n|
    threads << Thread.new(n) do |n|
      queue.pop
      system "convert -layers optimize -loop 0 -delay 10 -dispose 2 -transparent white plot/#{n}/*.png plot/dv_#{n}.gif"     if input[2]
      system "convert -layers optimize -loop 0 -delay 10 -dispose 2 -transparent white plot/o#{n}/*.png plot/orbit_#{n}.gif" if not input[2]
      queue.push 0
    end
  end
  threads.each{|t| t.join}

when "v"
  exit unless input[1]
  dtx  = 0.001
  dty  = 0.001
  vp   = [0, 0, 0]

  threads = []
  queue   = Queue.new
  10.times{|i| queue.push 0}

  idata = File.read "plot_dv.plt"
  index = input[1].to_i
  FileUtils.mkdir_p "plot/#{index}"
  [*0..99].each do |n|
    # next if n % 10 > 0
    sleep 0.2
    threads << Thread.new(n) do |n|
      queue.pop
      p n if n % 10 == 0
      i = 200
      j = n * 2

      i = n
      j = (i * 0.54 + 38).round
      # i = n
      # j = (i * (65 - 55) / 100.0 + 55).round

      of   = "plot/#{index}/dv_%03d.png" % n
      data = setterm  idata, of
      s    = "file_data using 1:2:($3<30?$3:30) every ::#{j}:#{i}:#{j}:#{i} w p pt 7 ps %.1f lc rgb '%s'"
      data = setmanual data, "replot %s, %s\n" % [s % [3, "black"], s % [2, "white"]]
      data.sub! /set pm3d.*$/, "set pm3d map" if index % 2 == 0
      data.sub! /set pm3d.*$/, "set pm3d"     if index % 2 == 1
      # data.sub! (/^(replot.*)$/){"# #{$1}"}
      plot data, n
      # system "convert -resize 400 #{of} #{of}"
      queue.push 0
    end
  end
  threads.each{|t| t.join}

when "orbit"
  system "start wgnuplot plot_orbit.rb"

when "edit"
  system "start plot_result.txt"

when "pt"
  col = {0=>"red", 1=>"orange", 2=>"green", 3=>"purple", 4=>"cyan"}
  100.times do |i|
    data = File.read "plot_debris.rb"
    tex  = (0..3).map{|j| "replot file_out2 every :::#{i*4+j}::#{i*4+j} with lines linewidth 2.0 linecolor rgb '#{col[j]}' notitle"}.join("\n")
    # p tex
    # exit
    data.gsub! /^# x$/, tex
    data.gsub! "points1.txt", "points#{i + 1}.txt"
    File.write "plot_debris2.rb", data
    system "wgnuplot plot_debris2.rb"
    # exit
    File.rename "orbit.png", "3/orbit%03d.png" % i
  end
end

time_end = Time.now

if time_end - time_start > 10
  puts "", "-" * 80
  puts "ARGS   : #{$*}"
  puts "START  : %s" % time_start.strftime("%Y/%m/%d %T")
  puts "FINISH : %s" % time_end.strftime("%Y/%m/%d %T")
  puts "TIME   : %.2f s" % (time_end - time_start)
end

__END__

行刻み
ブロック刻み
初期行
初期ブロック
終了行
終了ブロック
