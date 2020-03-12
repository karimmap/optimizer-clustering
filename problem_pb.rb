# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: problem.proto

require 'google/protobuf'

Google::Protobuf::DescriptorPool.generated_pool.build do
  add_file("problem.proto", :syntax => :proto3) do
    add_message "problem.Matrix" do
      repeated :time, :float, 1
    end
    add_message "problem.Location" do
      optional :longitude, :float, 1
      optional :latitude, :float, 2
      optional :matrix_index, :uint32, 3
    end
    add_message "problem.Service" do
      repeated :quantities, :float, 1
      optional :duration, :uint32, 2
      optional :matrix_index, :uint32, 3
    end
    add_message "problem.Capacity" do
      optional :limit, :float, 1
    end
    add_message "problem.Vehicle" do
      optional :id, :string, 1
      repeated :capacities, :message, 2, "problem.Capacity"
      optional :start_location, :message, 4, "problem.Location"
      optional :end_location, :message, 5, "problem.Location"
      optional :duration, :float, 6
      repeated :day_indices, :int32, 7
      optional :matrix_id, :string, 8
    end
    add_message "problem.Problem" do
      repeated :vehicles, :message, 1, "problem.Vehicle"
      repeated :services, :message, 2, "problem.Service"
      repeated :matrices, :message, 3, "problem.Matrix"
    end
  end
end

module Problem
  Matrix = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Matrix").msgclass
  Location = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Location").msgclass
  Service = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Service").msgclass
  Capacity = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Capacity").msgclass
  Vehicle = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Vehicle").msgclass
  Problem = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("problem.Problem").msgclass
end
