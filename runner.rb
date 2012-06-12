require "awesome_print"
require "vienna_rna"

sequence              = ARGV.first || (raise ArgumentError.new("No RNA sequence provided"))
generate_command_only = ARGV.last =~ /command/i

puts "NOTE: If scaling is not set to 1, solutions will need to be unscaled in Matlab."
puts "Making RNAbor..."
puts %x|make|
puts "Running RNAbor with #{sequence}..."
puts (output = %x|./RNAbor -s #{sequence}|)
puts "Parsing output..."

points = ViennaRna::Rnabor.parse(output) { |line| eval("Complex(#{line})") }

puts "Parsed roots and solutions:"
ap points

matrix  = "[%s]" % points.map(&:first).map { |root| (0...points.length).map { |power| root ** power }.join(" ") }.join("; ")
vector  = "[%s]" % points.map(&:last).join("; ")
command = case ARGV.last
when /.*fft.*/ then
  "b = %s; z = fft(b); z1 = z / real(z(1))" % vector
else
  "A = %s; b = %s; x = A \\ b" % [matrix, vector]
end

if generate_command_only
  puts "The Matlab command is:"
  
  case ARGV.last
  when /.*fft.*/ then
    puts "b = %s;" % vector
  else
    puts "A = %s;\n\n\n" % matrix
    puts "b = %s;\n\n\n" % vector
    puts command
  end
else
  puts "Going to run the following command in Matlab:"
  puts command

  %x|matlab -r "#{command}"|
end