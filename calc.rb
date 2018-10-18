#! ruby -Ku
require 'fileutils'
require 'socket'

def my_address
 udp = UDPSocket.new
 udp.connect "128.0.0.0", 7
 adrs = Socket.unpack_sockaddr_in(udp.getsockname)[1]
 udp.close
 adrs
end

DIR_IO   = "../dat"
DIR_WORK = "/path/to/worker"
TABLE    = ["192.168.10.10", "192.168.10.11", "192.168.10.13"]
# TABLE    = ["133.10.98.155"]
HOST     = my_address# and "yunibo.dip.jp"
ADDR     = 36700
PARA     = [30, 20, 0]
PARA_A   = PARA.map.with_index{|e, i| (0...e).map{|n| [i, n]}}.inject(:+).select{|n| n[0] < 2}

case HOST
when TABLE[0]
  NUM_PARA = PARA[0]
when TABLE[1]
  NUM_PARA = PARA[1]
when TABLE[2]
  NUM_PARA = PARA[2]
else
  NUM_PARA = PARA[0]
end

def calc input
  File.write "../dat/ADR_GA_Para.in",  input[0]
  File.write "../dat/ELEMENTS_JPL.in", input[1]

  system "ADR_GA.exe > NUL"

  File.read("../dat/ARDR_GA.out").scan(/^\s*\d.+$/)[0].split(/\s+/)[-3, 2].join(" ")
end

def open_server n
  host = HOST
  port = ADDR + n
  sleep n * 0.1
  puts "server start #{n} #{Time.now}"

  TCPServer.open(host, port) do |server|
    loop do
      client = server.accept
      input  = client.gets
      break if input =~ /close/
      result = calc input.chomp.split "#"
      client.puts result
      client.close
    end
  end
  sleep n * 0.1
  puts "server close #{n} #{Time.now}"
end

def start_server pn = NUM_PARA
  threads = []
  pn.times do |i|
    bin = "#{DIR_WORK}/#{i}/bin"
    dat = "#{DIR_WORK}/#{i}/dat"
    exe = "#{bin}/ADR_GA.exe"

    FileUtils.mkdir_p bin
    FileUtils.mkdir_p dat
    FileUtils.cp_r("ADR_GA.exe", exe) unless File.exist? exe

    threads << Thread.new(i) do |j|
      system "ruby #{__FILE__} open_server #{j}"
    end
  end
  threads.each{|t| t.join}
end

def close
  threads = []
  PARA_A.each do |e|
    threads << Thread.new(e) do |a|
      begin
        submit *a, "close"
      rescue
        puts "ERROR: #{a}"
      end
    end
  end
  threads.each{|t| t.join}
end

def submit dn, pn, input
  addr   = TABLE[dn]
  port   = ADDR + pn
  # p [addr, port]
  socket = TCPSocket.open addr, port
  socket.puts input
  result = socket.gets
  socket.close
  result and result.chomp
end

def main
  threads = []
  locks   = Queue.new
  qout    = Queue.new

  PARA_A.each{|e| locks.push e}
  qout.push []

  order = (File.readlines "#{DIR_IO}/order.txt").map(&:chomp)
  all   = order.size

  all.times do |i|
    threads << Thread.new(i, order[i]) do |j, o|
      begin
        e   = locks.pop
        res = submit *e, o
        # res = `ruby #{__FILE__} run #{n} \"#{o}\"`
        locks.push e

        out    = qout.pop
        out[j] = res
        qout.push out
      rescue
        puts "ERROR: #{j}"
      end
    end
  end

  threads.each{|t| t.join}

  File.open("#{DIR_IO}/res.txt", "w") do |file|
    file.puts qout.pop
  end
  all
end

case $*[0]
when "run"
  main

when "test"
  total = 0
  start = Time.now
  20.times do
    p total += main
  end
  finish = Time.now
  elap = finish - start
  p [elap, total / elap]

when "start"
  start_server

when "close"
  close

when "open_server"
  Dir.chdir("#{DIR_WORK}/#{$*[1]}/bin") do
    open_server $*[1].to_i
  end
  # e = 3
  # start = Time.now
  # (10**e).times do |i|
  #   res = `echo #{i} #{i+1} #{i+2} | a`
    # res = `a.exe #{i} #{i+1} #{i+2}`
    # p res.scan(/\d+/).map(&:to_i).inject(:+)
  # end
  # finish = Time.now
  # elap = finish - start
  # p [elap, 10**e / elap]

when "close0"
  PARA_A[0..19].each do |e|
    submit *e, "close"
  end

when "run0"
  Dir.chdir("#{DIR_WORK}/#{$*[1]}/bin") do
    print calc $*[2].split "#"
  end
end

__END__
