#! ruby -Ku

path = "C:/yam/Studies/Programs/Fortran/FEC/src/lib/*.f08"
path = "../src3/*.f08"

puts "total: %i" % Dir[path].map { |file|
  File.readlines(file).size.tap{ |s| puts "%-20s: %i" % [File.basename(file, ".*"), s] }
}.sum


__END__
