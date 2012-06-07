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
      table.decimal :scaling_factor
      table.string  :algorithm
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
      table.decimal    :real,      precision: 63, scale: 15
      table.decimal    :imaginary, precision: 63, scale: 15
    end  
  end
end

unless ActiveRecord::Base.connection.execute("show tables").map(&:this).flatten.include?("numeros")
  BuildNumeros.up
end

class Run < ActiveRecord::Base
  has_many :arguments, class_name: "Numero"
  has_many :solutions, class_name: "Numero"
  
  validates_presence_of :sequence, :scaling_factor, :algorithm
end

class Numero < ActiveRecord::Base
  belongs_to :run
end

class RealNumero < Numero
  validates_presence_of :real
  
  def value
    real
  end
end

class ComplexNumero < Numero
  validates_presence_of :real, :imaginary
  
  def value
    Complex(real, imaginary)
  end
end