require "awesome_print"

sequence              = ARGV.first || (raise ArgumentError.new("No RNA sequence provided"))
generate_command_only = ARGV.last =~ /command/i

puts "NOTE: If scaling is not set to 1, solutions will need to be unscaled in Matlab."
puts "Making RNAbor..."
puts %x|make|
puts "Running RNAbor with #{sequence}..."
puts (output = %x|./RNAbor -s #{sequence}|)
puts "Parsing output..."

points = output.split(/\n/).reject { |line| 
  line.empty?
}.drop_while { |line|
  line !~ /^roots and solutions/i
}.reverse.drop_while { |line|
  line !~ /solution/i
}.reverse[1..-2].map { |line|
  line.strip.split(/\s\s+/).map do |complex|
    eval("Complex(#{complex})")
  end
}

puts "Parsed roots and solutions:"
ap points

matrix  = "[%s]" % points.map(&:first).map { |root| (0...points.length).map { |power| root ** power }.join(" ") }.join("; ")
vector  = "[%s]" % points.map(&:last).join("; ")
command = "A = %s; b = %s; x = A \\ b" % [matrix, vector]

if generate_command_only
  puts "The Matlab command is:"
  puts command
else
  puts "Going to run the following command in Matlab:"
  puts command

  %x|matlab -r "#{command}"|
end