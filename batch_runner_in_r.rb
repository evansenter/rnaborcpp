require "./driver.rb"
require "vienna_rna"
require "awesome_print"

Run.where(algorithm: :complex_me).each do |template_run|
  ap (run = Run.create(sequence: template_run.sequence, scaling_factor: template_run.scaling_factor, algorithm: ARGV.last))
  
  template_run.arguments.each do |template_argument|
    run.arguments << ComplexArgument.new(k: template_argument.k, value: template_argument.value)
  end
  
  template_run.solutions.each do |template_solution|
    run.solutions << ComplexSolution.new(k: template_solution.k, value: template_solution.value)
  end
  
  vector   = "c(%s)" % run.solutions.map(&:pretty).join(", ")
  command  = "Rscript -e 'a <- #{vector}; b <- fft(a) / Re(fft(a)[1]); lapply(b, function(i) { sprintf(\"%g, %g\", Re(i), Im(i)) })'"
  response = %x|#{command}|
  
  response.split(/\n/).reject(&:empty?).select do |line| 
    line =~ /^\[1\] "(.+)"$/
  end.map do |line| 
    line.match(/^\[1\] "(.+)"$/)[1] 
  end.map do |number| 
    eval("Complex(#{number})")
  end.each_with_index.map do |count, k|
    run.counts << ComplexCount.new(k: k, value: count)
  end
end
