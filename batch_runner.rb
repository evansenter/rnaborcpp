require "./driver.rb"
require "vienna_rna"
require "awesome_print"

template = "GGGGGCCCCC" * 10

(10.step(50, 10).to_a + 55.step(75, 5).to_a + [100]).map do |sequence_length|
  (1.0.step(2.0, 0.5).to_a + 2.25.step(3.0, 0.25).to_a).map do |scaling_factor|
    sequence = template[0, sequence_length]
    
    ap (run = Run.create(sequence: sequence, scaling_factor: scaling_factor, algorithm: ARGV.last))
    
    results = ViennaRna::Rnabor.new(sequence).run(scaling_factor: scaling_factor)
    
    run.update_attributes(count: results.parse_total_count)
    
    results.parse_points.each_with_index do |point_array, i|
      run.arguments << ComplexArgument.new(k: i, value: point_array.first)
      run.solutions << ComplexSolution.new(k: i, value: point_array.last)
    end
    
    results.parse_counts.each do |k, count|
      run.counts << RealCount.new(k: k, value: count)
    end
  end
end