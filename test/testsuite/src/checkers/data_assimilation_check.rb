#!/usr/bin/env ruby

# COSMO TECHNICAL TESTSUITE
#
# This script checks whether the run assimilated the correct
# number of observations and information on the number of input
# observations (cdfin files) and output observations (feedobs files)

# Author       Daniel Leuenberger
# Email        cosmo-wg6@cosmo.org
# Maintainer   daniel.leuenberger@meteoswiss.ch

require "pp"
require "#{File.expand_path(File.dirname(__FILE__))}/yu.rb"

def exit20_on_condition(condition,message,verbose=2)
  if condition
    if verbose > 0
      puts "Configuration:"	
      puts "TS_VERBOSE   = #{ENV["TS_VERBOSE"]}"
      puts "TS_RUNDIR    = #{ENV["TS_RUNDIR"]}"
      puts "TS_REFOUTDIR = #{ENV["TS_REFOUTDIR"]}"
      puts "TS_LOGFILE   = #{ENV["TS_LOGFILE"]}"
      puts message
    end 
    exit(20)
  end
end

def check_env
  verbose   = ENV["TS_VERBOSE"]
  exit20_on_condition(verbose.nil?,"TS_VERBOSE is undefined")
  rundir    = ENV["TS_RUNDIR"]
  exit20_on_condition(rundir.nil? || !File.directory?(rundir),"TS_RUNDIR is undefined or not a directory")
  refoutdir = ENV["TS_REFOUTDIR"]
  exit20_on_condition(refoutdir.nil? || !File.directory?(refoutdir),"TS_REFOUTDIR is undefined or not a directory")
  stdout    = ENV["TS_LOGFILE"]
  exit20_on_condition(stdout.nil?,"TS_LOGFILE is undefined")
  return {:verbose=>verbose.to_i,:rundir=>rundir.to_s,:refoutdir=>refoutdir.to_s,:stdout=>stdout.to_s}
end

def process_stdout(file,verbose)

  lines = File.open(file,"r").readlines.map{|l|l.chomp}
  out = {}


  #### remove this check once single precision is fixed for data assimilation !
  puts "check for single precision execution" if verbose > 1
  pattern = /RUNNING IN SINGLE PRECISION/
  unless lines.select{|x|x=~pattern}.empty?		
    puts "currently, running in single precision is not supported in data assimilation, skip" if verbose > 0
    exit(15)
  end 
  #### end remove
 
  puts "check number of time steps" if verbose > 1
  pattern = /STEP *\d+ *\d*/
  out[:number_of_timesteps] = lines.select{|x|x=~pattern}.length
  
  puts "check obs read from cdfin" if verbose > 1
  pattern = /processing *(.*)from file cdfin(.*)/
  out[:cdfin_obs] = lines.select{|x|x=~pattern}.map{|x|
    x=~pattern
    "#{$1}#{$2}"
  }

  puts "check information on feedobs file" if verbose > 1
  pattern = /feedobs file: *(.*),(.*),(.*),(.*).*$/
  out[:feedobs_reports] = lines.select{|x| x=~pattern}.map{|x|
    x=~pattern
    "#{$1}#{$2}#{$3}#{$4}"
  }
  puts "check information on assimilated reports" if verbose > 1
  pattern = /REPORTS WITH OBS INCREMENTS:\s*(\d+)/
  out[:reports_assimilated] = lines.select{|x|x=~pattern}.map{|x|
    x=~pattern
    "REPORTS:#{$1}"
  }
  pattern = /(\d+)\s*multi-level:\s*(\d+)\s*TEMP,\s*(\d+)\s*PILOT,\s*(\d+)\s*WINDPROF,\s*(\d+)\s*RADAR-VAD/
  out[:reports_assimilated] << lines.select{|x|x=~pattern}.map{|x|
    x=~pattern
    "multi-level:#{$1},TEMP:#{$2},PILOT:#{$3},WINDPROF:#{$4},RADAR-VAD:#{5}"
  }
  pattern = /\s*(\d+)\s*AIRCRAFT,\s*(\d+)\s*RASS,\s*(\d+)\s*GPS,\s*(\d+)\s*RETRIEVALS/
  out[:reports_assimilated] << lines.select{|x|x=~pattern}.map{|x|
    x=~pattern
    "AIRCRAFT:#{$1},RASS:#{$2},GPS:#{$3},RETRIEVALS:#{$4}"
  }
  
  puts "check information on analysis increments" if verbose > 1
  pattern = /AI-box: *(.*)$/
  out[:analysis_increments] = lines.select{|x|x=~pattern}.map{|x|
    x=~pattern
    $1
  }
  return out
end


def process_yustats(file,verbose)
  y = Yustats::YuStats.new(file)
  out = {}
  y.data.each_pair do |ot,a|
    next if ot.nil?
    out[ot.to_sym]="observations processed:#{a[0]}, active:#{a[1]},passive:#{a[2]},rejected:#{a[3]}"
  end
  return out
end


def compare(h_ref,h_exp,verbose)
  status = ""
  h_ref.each_pair do |k,v|
    if v != h_exp[k]
      if verbose > 0
        puts "info about #{k.to_s} differs"
        if verbose > 1
          puts "reference:"
          puts v
          puts "experiment:"
          puts h_exp[k]
        end
      end
      status = "fail"
    end
  end
  exit20_on_condition(status == "fail","",verbose)
end

###########  Main script starts here  #####################################

# test environment variables
#ENV["TS_VERBOSE"] = "2"
#ENV["TS_RUNDIR"] = Dir.pwd+"/exp/"
#ENV["TS_REFOUTDIR"] = Dir.pwd+"/ref"
#ENV["TS_LOGFILE"] = "logfile"

# check environment and set verbose status
env = check_env()
verbose = env[:verbose]

#### stdout ####

puts "check stdout..." if verbose > 1
stdout_ref = env[:refoutdir]+"/"+env[:stdout]
stdout_exp = env[:rundir]   +"/"+env[:stdout]
puts "process reference  file #{stdout_ref}" if verbose > 1
puts "process experiment file #{stdout_exp}" if verbose > 1
out_ref = process_stdout(stdout_ref,verbose)
out_exp = process_stdout(stdout_exp,verbose)

# skip DA checker, if number of time steps differ between reference and experiment
# since the number of assimilated observations depends on the number of time
# steps
if (out_ref[:number_of_timesteps] != out_exp[:number_of_timesteps])
  puts "number of time steps differ, skip" if verbose > 0
  exit(15)
end

# compare information of stdout
compare(out_ref,out_exp,verbose)

#### YUSTATS ####

puts "check YUSTATS..." if verbose > 1
yustats_ref = env[:refoutdir]+"/YUSTATS"
yustats_exp = env[:rundir]   +"/YUSTATS"
puts "process reference  file #{yustats_ref}" if verbose > 1
puts "process experiment file #{yustats_exp}" if verbose > 1
out_ref = process_yustats(yustats_ref,verbose)
out_exp = process_yustats(yustats_exp,verbose)
compare(out_ref,out_exp,verbose)

# all went well, so yippieeeee and exit with 0
exit(0)
