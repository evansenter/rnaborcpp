require "./driver.rb"
require "vienna_rna"
require "awesome_print"

module Enumerable
  def sum
    inject { |sum, i| sum + i }
  end
end

sequence       = ARGV[0]
scaling_factor = ARGV[1].to_f
precision      = ARGV[2].to_f
results        = ViennaRna::Rnabor.new(sequence).run(scaling_factor: scaling_factor)
inferred_c     = results.parse_total_count ** (1.0 / (sequence.length - 1))

puts "RUNNING WITH NO SCALING IN RECURSIONS\n\n\n"
puts results.response

vector   = "c(%s)" % results.parse_points.map(&:last).map { |y| 10 ** precision * y / results.parse_total_count }.join(", ")
command  = "Rscript -e 'a <- #{vector}; fft(a) / length(a);'"
response = %x|#{command}|

puts command
puts response

# results = response.split(/\n/).map { |line| line.strip.match(/\[\d+\]\s+(.*)$/)[1].split(/\s+/) }.flatten.map { |i| i.match(/(-?\d+\.\d+e[\+-]\d+)/)[1].to_f.truncate / 10.0 ** precision }
results = response.split(/\n/).map { |line| line.strip.match(/\[\d+\]\s+(.*)$/)[1].split(/\s+/) }.flatten.map { |i| i.match(/(-?\d+\.\d+e[\+-]\d+)/)[1].to_f * results.parse_total_count }

ap results
ap results.sum

puts "RUNNING WITH SCALING IN RECURSIONS\n\n\n"
puts "Inferred scaling factor: #{inferred_c}"

results = ViennaRna::Rnabor.new(sequence).run(scaling_factor: inferred_c)

puts results.response

vector   = "c(%s)" % results.parse_points.map(&:last).map { |y| 10 ** precision * y }.join(", ")
command  = "Rscript -e 'a <- #{vector}; fft(a) / length(a);'"
response = %x|#{command}|

puts command
puts response

# results = response.split(/\n/).map { |line| line.strip.match(/\[\d+\]\s+(.*)$/)[1].split(/\s+/) }.flatten.map { |i| i.match(/(-?\d+\.\d+e[\+-]\d+)/)[1].to_f.truncate / 10.0 ** precision }
results = response.split(/\n/).map { |line| line.strip.match(/\[\d+\]\s+(.*)$/)[1].split(/\s+/) }.flatten.map { |i| i.match(/(-?\d+\.\d+e[\+-]\d+)/)[1].to_f * results.parse_total_count }

ap results
ap results.sum