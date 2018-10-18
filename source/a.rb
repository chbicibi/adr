#! ruby -Ku
st = Time.now
system "a.exe"
et = Time.now
puts "%.2f s" % (et - st)

__END__
