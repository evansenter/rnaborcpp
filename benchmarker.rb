require "rubygems"
require "mysql2"
require "active_record"
require "vienna_rna"
require "awesome_print"

class Object; def this; self; end; end

ActiveRecord::Base.establish_connection(config = { adapter: "mysql2", username: "root", reconnect: true })

unless ActiveRecord::Base.connection.execute("show databases").map(&:this).flatten.include?("fftbor_performance")
  ActiveRecord::Base.connection.create_database("fftbor_performance")
end

ActiveRecord::Base.establish_connection(config.merge(database: "fftbor_performance"))

class BuildRun < ActiveRecord::Migration
  def self.up
    create_table :runs do |table|
      table.string  :sequence
      table.integer :sequence_length
      table.string  :algorithm
      table.decimal :time, precision: 20, scale: 3
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

class Run < ActiveRecord::Base
  validates_presence_of :sequence, :sequence_length, :algorithm
  
  def self.generate_sequence(sequence_length)
    sequence_length.times.inject("") { |string, _| string + %w[A U C G][rand(4)] }
  end
end

ViennaRna.debug = false

20.step(260, 20).inject({}) { |hash, size| hash.merge(size => (case size; when 1..150 then 100; when 151..200 then 5; else 3; end)) }.each do |size, iterations|
  iterations.times do
    sequence = Run.generate_sequence(size)

    ap (fftbor = Run.create(sequence: sequence, sequence_length: size, algorithm: "fftbor", time: ViennaRna::Fftbor.new(sequence: sequence).run.runtime.real)).attributes
    ap (rnabor = Run.create(sequence: sequence, sequence_length: size, algorithm: "rnabor", time: ViennaRna::Rnabor.new(sequence: sequence).run.runtime.real)).attributes
  end
end