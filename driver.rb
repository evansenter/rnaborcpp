# In the future this should use double columns in MySQL for better range.

class Object; def this; self; end; end

require "rubygems"
require "mysql2"
require "active_record"

config = YAML::load(File.open("database.yml"))
ActiveRecord::Base.establish_connection(config)

unless ActiveRecord::Base.connection.execute("show databases").map(&:this).flatten.include?("rnabor")
  ActiveRecord::Base.connection.create_database("rnabor")
end

ActiveRecord::Base.establish_connection(config.merge(database: "rnabor"))

class BuildRun < ActiveRecord::Migration
  def self.up
    create_table :runs do |table|
      table.string  :sequence
      table.integer :sequence_length
      table.decimal :scaling_factor, precision: 2, scale: 1
      table.string  :algorithm
      table.float   :count
      table.timestamps
    end 
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("runs")
  BuildRun.up
end

class BuildNumeros < ActiveRecord::Migration
  def self.up
    create_table :numeros do |table|
      table.belongs_to :run
      table.string     :type
      table.integer    :k
      table.decimal    :real,      precision: 63, scale: 15
      table.decimal    :imaginary, precision: 63, scale: 15
      table.timestamps
    end  
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("numeros")
  BuildNumeros.up
end

class Run < ActiveRecord::Base
  has_many :arguments, class_name: "Numero", dependent: :destroy, conditions: "type like '%Argument'"
  has_many :solutions, class_name: "Numero", dependent: :destroy, conditions: "type like '%Solution'"
  has_many :counts,    class_name: "Numero", dependent: :destroy, conditions: "type like '%Count'"
  
  validates_presence_of :sequence, :sequence_length, :scaling_factor, :algorithm
  
  def total_estimated_error
    counts.inject(0.0) do |sum, count|
      value = (count.value.class == Complex ? count.value.real : count.value)
      
      sum + (value < 0 ? -value : (value.round - value).abs)
    end
  end
  
  def avg_estimated_error
    total_estimated_error / (sequence.length + 1)
  end
end

class Numero < ActiveRecord::Base
  belongs_to :run
  
  validates_presence_of :run_id
  
  def pretty
    value.to_s
  end
end

class RealNumero < Numero
  validates_presence_of :real
  
  def value
    real
  end
  
  def value=(real)
    self.real = real
  end
end

class ComplexNumero < Numero
  validates_presence_of :real, :imaginary
  
  def value
    Complex(real, imaginary)
  end
  
  def value=(complex)
    self.real      = complex.real
    self.imaginary = complex.imaginary
  end
end

class RealArgument < RealNumero; end
class RealSolution < RealNumero; end
class RealCount < RealNumero; end
class ComplexArgument < ComplexNumero; end
class ComplexSolution < ComplexNumero; end
class ComplexCount < ComplexNumero; end