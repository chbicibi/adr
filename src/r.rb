#! ruby -Ku
require 'fileutils'

EXE      = "adr.exe"
DIR_TASK = "task"
N        = 50

def start title = "none"
  return nil if Dir.exist? "#{DIR_TASK}/#{title}"

  (0...N).each do |i|
    return nil if File.file? "stop"
    dir = "#{DIR_TASK}/#{title}/#{i}"

    if i == 0
      FileUtils.mkdir_p dir
      FileUtils.cp EXE, "#{dir}/#{EXE}"
      FileUtils.cp "input.txt", "#{dir}/input.txt"
      FileUtils.cp_r "dat", "#{dir}/dat"
    else
      FileUtils.rm_r dir if Dir.exist? dir
      FileUtils.cp_r "#{DIR_TASK}/#{title}/0", dir
    end
  end

  (0...N).each do |i|
    dir = "#{DIR_TASK}/#{title}/#{i}"

    if (i % 10 < 9)
      system "start \"#{i}\" cmd /k ruby #{__FILE__} run #{dir} #{rand 2 ** 30}"
    else
      Dir.chdir(dir) do
        FileUtils.rm_f "stop"
        system "#{EXE} #{rand 2 ** 30}"
      end
    end
  end
end

def result
  return r2
  # dir = Dir["task/{p0,ch,test28}*"]
  dir = Dir["/path/to/result/{p0,ch,test28}*"]
  dir.each do |d|
    dn   = File.basename d
    list = Dir["#{d}/*/result.csv"]
    data = list.map do |f|
      (File.read f).scan(/^\d.+$/).map{|e| e.split(",")[1].strip.to_f}.uniq.sort[0]#.select{|e| e.to_i < 76}
    end.sort#.reverse
    best = data[0]
    mean = data.inject(:+) / data.size
    medi = data[data.size / 2]
    wost = data[-1]
    devi = (data.map{|e| (e - mean) ** 2}.inject(:+) / data.size) ** 0.5
    puts "%s: %2i B [%.5f] A [%.5f] M [%.5f] W [%.5f] D [%.5f]" % [dn, data.size, best, mean, medi, wost, devi]
    # puts "%.5f %.5f %.5f %.5f %.5f" % [best, mean, medi, wost, devi]
    # puts "#{dn}: #{best}"
  end
end

def r2
  dir = "ch2" # 1段階
  # dir = "p05" # 2段階
  # list = Dir["C:/Yamada/Society/20161210_進化計算シンポジウム/keep/#{dir}/*/result.csv"]
  list = Dir["/path/to/data/y/*/result.csv"]
  data = list.map do |f|
    n = f.scan(/\/(\d+)\//)[0][0].to_i
    [(File.read f).scan(/^\d.+$/).map{|e| e.split(",")[1].strip.to_f}.uniq.sort[0], n]#.select{|e| e.to_i < 76}
  end.sort#.reverse
  p *data
end

def search
  dir = Dir["src/*.f08"]
  kw  = /THash/
  rw  = "Hash"
  dir.each do |f|
    txt = File.read f
    hit = txt.scan(kw)
    next if hit.empty?
    p [f, *hit]
    # next
    File.write f, (txt.gsub kw, rw)
  end
end

# ------------------------------------------------------------------------------

input      = $*.reject{|e| e =~ /^\//}
opts       = $* - input
time_start = Time.now

if opts.include? "/c"
  FileUtils.rm_f Dir["src/*.{o,mod}"]
  system "cls"
end
if opts.include? "/m"
  Dir.chdir("src") do
    # FileUtils.rm_f Dir["*.o", "*.mod"]
    exit unless (p system "make" + (input[0] != "run" ? " #{input[0]}" : ""))
  end
end

case input[0]
when "start"
  start input[1]

when "run"
  if input[1].nil?
    FileUtils.rm_f "stop"
    system "#{EXE} #{rand 2 ** 30}"
  else
    (puts "nodir") or exit unless Dir.exist? input[1]
    Dir.chdir(input[1]) do
      FileUtils.rm_f "stop"
      system "#{EXE} #{input[2] or rand 2 ** 30}"
    end
  end

when "stop"
  (0...N).each do |i|
    File.write "#{DIR_TASK}/#{input[1]}/#{i}/stop", ""
  end

when "result"
  result

when "plot"
  system "start wgnuplot plot_pareto.rb"

when "pareto"
  system "start wgnuplot plot_pareto.rb"

when "orbit"
  system "start wgnuplot plot_orbit.rb"

when "edit"
  system "start plot_result.txt"

when "pt"
  col = {0=>"red", 1=>"orange", 2=>"green", 3=>"purple", 4=>"cyan"}
  100.times do |i|
    data = File.read "plot_debris.rb"
    tex  = (0..3).map{|j| "replot file_out2 every :::#{i*4+j}::#{i*4+j} with lines lw 2.0 lc rgb '#{col[j]}' notitle"}.join("\n")
    # p tex
    # exit
    data.gsub! /^# x$/, tex
    data.gsub! "points1.txt", "points#{i + 1}.txt"
    File.write "plot_debris2.rb", data
    system "wgnuplot plot_debris2.rb"
    # exit
    File.rename "orbit.png", "3/orbit%03d.png" % i
  end

when "search"
  search

when "test"
  system "test.exe #{rand 2 ** 30}"

else
  puts "wrong command"
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
