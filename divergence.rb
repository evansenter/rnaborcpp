require "rubygems"
require "mysql2"
require "active_record"
require "awesome_print"
require "vienna_rna"
require "diverge"

class Object; def this; self; end; end

ActiveRecord::Base.establish_connection(config = { adapter: "mysql2", username: "root", reconnect: true })

unless ActiveRecord::Base.connection.execute("show databases").map(&:this).flatten.include?("fftbor_divergence")
  ActiveRecord::Base.connection.create_database("fftbor_divergence")
end

ActiveRecord::Base.establish_connection(config.merge(database: "fftbor_divergence"))

class BuildRun < ActiveRecord::Migration
  def self.up
    create_table :runs do |table|
      table.string  :sequence
      table.integer :sequence_length
      table.string  :structure
      table.decimal :js_divergence, precision: 20, scale: 5
      table.decimal :fftbor_time,   precision: 20, scale: 3
      table.decimal :rnabor_time,   precision: 20, scale: 3
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

class Run < ActiveRecord::Base
  validates_presence_of :sequence, :sequence_length, :structure, :js_divergence, :fftbor_time, :rnabor_time
  
  def self.generate_sequence(sequence_length)
    sequence_length.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
  end
end

40.step(200, 20).each do |size|
  sequence = Run.generate_sequence(size)
         
  ["." * sequence.length, ViennaRna::Fold.new(sequence).run.structure].each do |structure|
    fftbor = ViennaRna::Fftbor.new(sequence: sequence, structure: structure).run
    rnabor = ViennaRna::Rnabor.new(sequence: sequence, structure: structure).run
    
    puts (fftbor_distribution = fftbor.parse_distribution).inspect
    puts (rnabor_distribution = rnabor.parse_distribution.map { |i| (i * 10 ** 6).truncate / 10.0 ** 6 }).inspect
    
    Run.create({
      sequence:        sequence, 
      sequence_length: size, 
      structure:       structure, 
      js_divergence:   Diverge.new(fftbor_distribution, rnabor_distribution).js,
      fftbor_time:     fftbor.runtime.real,
      rnabor_time:     rnabor.runtime.real
    })
  end
end