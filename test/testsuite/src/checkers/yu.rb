require 'pp'
require 'zlib'

# This ruby source file collects all ruby code necessary for the treatment
# of diagnostic ASCII output (YU Files) of the COSMO Model. The module YU
# stores information about all YU files and for each YU File, a specific ruby module
# is provided

module YU
  
  CTYPS = {
    11  => "MANUAL LAND SYNOP",
    14  => "AUTOM. LAND SYNOP",
    21  => "MANUAL SHIP SYNOP",
    22  => "ABBR. SHIP SYNOP",
    23  => "REDUCED SHIP SYNOP",
    24  => "AUTOM. SHIP SYNOP",
    41  => "CODAR AIRCRAFT",
    140 => "METAR",
    141 => "AIREP AIRCRAFT",
    241 => "CONST LEV BALLOON",
    144 => "AMDAR AIRCRAFT",
    244 => "ACAR  AIRCRAFT",
    63  => "BATHY SPHERE",
    64  => "TESAC",
    165 => "DRIFTING BUOY",
    35  => "LAND RADIO SONDE",
    36  => "SHIP RADIO SONDE",
    135 => "DROP SONDE",
    39  => "LAND ROCKET SONDE",
    40  => "SHIP ROCKET SONDE",
    32  => "LAND PILOT",
    33  => "SHIP PILOT",
    38  => "MOBILE PILOT",
    132 => "WIND PROFILER EUR",
    133 => "SODAR/RASS EUR",
    136 => "PROFILER/RASS USA",
    137 => "RADAR VAD",
    71  => "MSG-1 SEVIRI",
    206 => "NOAA-15 ATOVS",
    207 => "NOAA-16 ATOVS",
    208 => "NOAA-17 ATOVS",
    209 => "NOAA-18 ATOVS"
  }
  
  OBS_TYPES = {
    "synop"   => [11,14],
    "ship"    => [21,22,23,24],
    "aircraft"=> [41,141,241,144,244],
    "buoy"    => [165],
    "temp"    => [35,36],
    "pilot"   => [32,33],
    "wprof"   => [132]
  }
end

module Yuprmass
  
  HEADER  = "***** YUPRMASS begins *****"
  TRAILER = "***** YUPRMASS ends *****\\n\\n"
  
  class YuPrmass
    
    attr_reader :data
    attr_accessor :start
    
    
    def initialize(filename)
      @data = {}
      @start = nil
      str = nil
      File.open(filename) { |f|
        begin
          gz = Zlib::GzipReader.new(f)
          str = gz.read
        rescue
          f.rewind
          str = f.read
        ensure
          f.close unless f.nil?
        end 
      }
      parse(str) unless str.nil?
    end
    
    def parse(str)
      
      strbeg = str.index(HEADER)
      strbeg = 0 if strbeg.nil? 
      
      strend = str.index(TRAILER)
      strend = -1 if strend == nil
      buffer = str[strbeg..strend].split("\n")
      
      buffer.each{|line|
        a = line.split
        next if a.empty?
        ntstep = a[0].strip
        if a.length==12 and ntstep=~/^\d+$/
          @data[ntstep] = {
            "dpsdt"=>a[2].to_f,
            "ps"   =>a[3].to_f,
            "dse"  =>a[4].to_f,
            "mse"  =>a[5].to_f,
            "ke"   =>a[6].to_f,
            "vhmax"=>a[7].to_f,
            "wmax" =>a[8].to_f,
            "wa850"=>a[9].to_f,
            "wa500"=>a[10].to_f,
            "wa300"=>a[11].to_f
          }
        end
      }      
    end # parse
    
    def empty?
      @data.empty?
    end


  end # class
  
  
  
end # module

module Yustats
  
  HEADER  = "0 ---DISTRIBUTION OF PROCESSED/ACTIVE/PASSIVE/REJECTED REPORTS FOR"
  TRAILER = "***** YUSTATS ends *****\\n\\n"
  STATES  = {
    "processed" => 0,
    "active"    => 1,
    "passive"   => 2,
    "rejected"  => 3
  }
  
  class YuStats
    
    attr_accessor :data
    attr_accessor :start
    
    
    def initialize(filename)
      @data = {}
      @start = nil
      str = nil
      if File.file?(filename)
        File.open(filename) { |f|
          begin
            gz = Zlib::GzipReader.new(f)
            str = gz.read
          rescue
            f.rewind
            str = f.read
          ensure
            f.close unless f.nil?
          end 
        }
        parse(str) unless str.nil?
      end
    end
    
    def parse(str)
      
      strbeg = str.index(HEADER)
      return if strbeg.nil? 
      
      strend = str.index(TRAILER)
      strend = -1 if strend == nil
      buffer = str[strbeg..strend].split("\n")
      
      
      buffer.each{|line|
        if line =~ /code type\s+(\d+).*?\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/
          code_type,processed,active,passive,rejected = $1,$2.to_i,$3.to_i,$4.to_i,$5.to_i
          obs_name = YU::CTYPS[code_type.to_i]
          @data[obs_name] = [processed,active,passive,rejected]
        elsif line =~ /observation type *\d+\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/
          obs_name,processed,active,passive,rejected = $1,$2.to_i,$3.to_i,$4.to_i,$5.to_i
          @data[obs_name] = [processed,active,passive,rejected]
        end
      }      
      
      
    end # parse
    
    def empty?
      @data.empty?
    end


  end # class
  
  class LogInteg < YuStats
    def parse(str)
      buffer = str.split("\n")
      buffer.each do |line|
        if line =~ /lm_dwh_retrieve\((.*)\)\]: Successfully retrieved (\d+) messages from DWH/
          obs_name,number=$1,$2
          @data[obs_name] = number
        end
      end
    end
  end
  
end # module



# YUVERIF VOF files are diagnostic output files from the COSMO model.
#
# The format of the yuverif files is difficult to use with computer systems
#
# This module provides Class and methods to parse, sort, reformat Yuverif files.
# Additional function that can be used to extract parts of the information found
# in the Yuverif file (find_all) are also available
#
# For more information about the YUVEIRF VOF output of the COSMO model
# the reader should have a look at Section 8: Model Output in Part VII - User's Guide
# of the COSMO Model documentation
#
#
# This modules defines the YuVerif class.
#
#	YuVerif ->
#				:info, :runs, :report
#
#   @info is a Hash with the text of the Yuverif Header part
# 	@run is a Hash with information about the model runs stored in the Yuverif
#	@report is an array containing the reports that were found in the yuverif file
#
# In order to store report information this module defines the Report class:
#
# A report consists of two, possibly three parts:
#
#	Report ->
#			  	:repHeader, :basetime, :runs, :multi_level, :single_level, \
#  				:report_type, :regRepBody, :optRepBody, :nlevels
#
# repHeader is of class ReportHeader, regRepBody is of class RegRepBody,
# optRepBody (if available) is a Hash of Class OptRepBody (one entry per run)
# 
# Sample session:
#
# >> require 'yuverif.rb'
# >> yuverif = Yuverif::YuVerif.new
# >> yuverif.push_file('MyYuverif.file')
# >> yuverif.lenght (returns the number of reports)
# >> yuverif.datatypes (returns a Hash with the data types found)
# >> yuverif[3] (returns report number 4)
# >> yuverif.find_all(|report| report.code_type == 137} (returns a new yuverif with only VAD reports (code_type = 137)
# >> yuverif.write_txt('myFileName') 	creates ASCII outputs of the Yuverif
#										(one file per data type; Named: myFileName-yuverif-???.data.gz)
# etc...
#
# Authors:
# 	Christophe Hug (huc) MeteoSwiss / MO 2008
# 	Daniel Leuenberger (led) MeteoSwiss / MO 2008
# 	Oliver Marchand (mao) MeteoSwiss / MO 2006 (Original design of the VerifData class)
# 
module Yuverif

  MISSING = -9999
  
  HEADER  = "VOF: Verification Observation File:"
  TRAILER = "\n  -9 xxxxx"
  NUMBEROFRUNS = "Number of model runs to compare with observations:"
  ModelRunTypes = "no\.\|type\| of the run  \| period [hrs] \|\[1\/deg\.\]\| levels \| the model run"
  
  REPORT_TYPE_EXPECTED_COLUMNS = {
    1 => ["multi-level report",10,5,3],
    0 => ["complete surface report",22,14,1],
    -1 => ["short surface report",16,9,2],
    -2 => ["very-short surface-level report",11,5,3],
    -3 => ["upper air single report",11,4,4],
    -4 => ["GPS report",11,1,12]
  }
  
  PRINTORDER = {
    0 => "id",
    1 => "lat",
    2 => "lon",
    3 => "reltime",
    4 => "time",
    5 => "sta_height",
    6 => "mod_height",
    7 => "obs_type",
    8 => "cod_type",
    9 => "type",
    10 => "status",
    11 => "report_flags",
    12 => "station_char",
    13=> "data"
  } 
  
  
  def Yuverif.strarrtonumarr(sarr, scale=1.0, missing=9999,offset=0.0)
    sarr.collect!{|elem| Yuverif.convert(elem,scale,missing,offset)}
  end
  
  def Yuverif.convert(number,scale=1.0,missing=9999,offset=0.0)
    if number.to_i == missing
      return MISSING
    else
      return number.to_f*scale
    end
  end
  
  ############################################################
  
  #
  # This class stores the content of a Yuverif file but can also
  # be use to store a collection of reports selected by the self.find_all
  # method
  #
  class YuVerif
    
    attr_reader :info, :runs, :reports
    
    
    def initialize
      @info = []
      @runs = []
      @reports = []
    end
    
    def push_file(filename)
      str = nil
      File.open(filename) { |f|
        begin
          gz = Zlib::GzipReader.new(f)
          str = gz.read
        rescue
          f.rewind
          str = f.read
        ensure
          f.close unless f.nil?
        end }
      parse(str) unless str.nil?
    end
    
    def write_bin(filename)
      file = File.open(filename,"w")
      file.write(Marshal.dump(self))
      file.close
    end
    
    def push_bin(filename)
      newYuverif = Marshal.load(File.read(filename))
      @info = newYuverif.info
      @runs = newYuverif.runs
      @reports = newYuverif.reports
    end
    
    def basetime
      @info["basetime"]
    end
    
    def number_of_runs
      @runs["start"].length
    end
    
    def length
      @reports.length
    end
    
    def info=(newval)
      @info = newval
    end
    
    def runs=(newval)
      @runs = newval
    end
    
    def push_report(newval)
      newval.inject(@reports) { |reports,elem|
        reports << if elem.class != Yuverif::Report then
                     puts "ERROR:: You can only add elements of Class Yuverif::Report"
                     nil
                   else
                     elem
                   end
      }
      @reports.compact! # remove nil elements
    end
    
    def [](elem)
      @reports[elem]
    end
    
    def find_all
      newYuverif = Yuverif::YuVerif.new
      newYuverif.info = self.info.dup
      newYuverif.runs = self.runs.dup
      newReports = self.reports.find_all {|myElem| yield(myElem)}
      newYuverif.push_report(newReports)
      return newYuverif
    end
    
    def write_txt(filename, format='old')
      
      ## We write the reports and the increments using the format:
      ## Number_of_runs basetime start end
      ## 1 initial end mesh nlevels description
      ## 2 initial end mesh nlevels description
      ## ... initial end mesh nlevels description
      ## n initial end mesh nlevels description
      ## StationID
      ## Report header (with time)
      ## Report Body (on one single line) 1
      ## Report Body (on one single line) 2
      ## Report Body (on one single line) ...
      ## Report Body (on one single line) n            
      ## Report increments (on one single line) 1
      ## Report increments (on one single line) 2
      ## Report increments (on one single line) ...
      ## Report increments (on one single line) n              
      ## 
      ## The reports are written sorted according to the station ID first and time.
      
      

      str_header = []
      n_of_runs = self.number_of_runs
      
      # Text file header
      str_header << ["#{"%3d" % n_of_runs} #{self.basetime.utc.strftime("%Y%m%d%H")} #{self.info["start"]} #{self.info["end"]}"]

      ([] <<
           (1..n_of_runs).collect{|val| "%3d" % val } <<
           @runs["start"].collect{|val| val.utc.strftime(" %Y%m%d%H")} <<
           @runs["endtime"].collect{|val| "%8.4f" % val } <<
           @runs["mesh_width"].collect{|val| "%6.2f" % val } <<
           @runs["nlevels"].collect{|val| " %3d" % val } <<
           @runs["type"].collect{|val| " %2d" % val } <<
           @runs["description"].collect{|val| " " + val } ).transpose.each {|val| str_header << [val.join]}
      str_header << [""] # line break
      
      # Reports
      str = {}
      @reports.sort_by{|reportkey| 
        sorting_key = reportkey.repHeader.station_id.dup.to_s
        sorting_key << ":" + reportkey.time.strftime("%Y%m%d%H%M")
      }.each{ |report|
        code_type = report.repHeader.code_type
        if str.has_key?(code_type) then
          str[code_type] << [report.write_txt(format)]
        else
          str[code_type] = []
          str[code_type] << [report.write_txt(format)]
        end }
      
      # Write everything to disc
      str.each{ |key, value|
        myFile = Zlib::GzipWriter.open("#{filename}.yuverif-#{"%03d" % key}.data.gz")
        myFile.write(str_header.join("\n"))
        myFile.write(value.join("\n"))
        myFile.close
      }
      
    end

   def write_out(filename)
      
      ## We write the reports and the increments using the format:
      ## Number_of_runs basetime start end
      ## 1 initial end mesh nlevels description
      ## 2 initial end mesh nlevels description
      ## ... initial end mesh nlevels description
      ## n initial end mesh nlevels description
      ## StationID
      ## Report header (with time)
      ## Report Body (on one single line) 1
      ## Report Body (on one single line) 2
      ## Report Body (on one single line) ...
      ## Report Body (on one single line) n            
      ## Report increments (on one single line) 1
      ## Report increments (on one single line) 2
      ## Report increments (on one single line) ...
      ## Report increments (on one single line) n              
      ## 
      ## The reports are written sorted according to the station ID first and time.
      
      

      str_header = []
      n_of_runs = self.number_of_runs
      
      # Text file header
      str_header << ["#{"%3d" % n_of_runs} #{self.basetime.utc.strftime("%Y%m%d%H")} #{self.info["start"]} #{self.info["end"]}"]

      ([] <<
           (1..n_of_runs).collect{|val| "%3d" % val } <<
           @runs["start"].collect{|val| val.utc.strftime(" %Y%m%d%H")} <<
           @runs["endtime"].collect{|val| "%8.4f" % val } <<
           @runs["mesh_width"].collect{|val| "%6.2f" % val } <<
           @runs["nlevels"].collect{|val| " %3d" % val } <<
           @runs["type"].collect{|val| " %2d" % val } <<
           @runs["description"].collect{|val| " " + val } ).transpose.each {|val| str_header << [val.join]}
      str_header << [""] # line break
      
      # Reports
      str = {}
      @reports.sort_by{|reportkey| 
        sorting_key = reportkey.repHeader.station_id.dup.to_s
        sorting_key << ":" + reportkey.time.strftime("%Y%m%d%H%M")
      }.each{ |report|
        code_type = report.repHeader.code_type
        if str.has_key?(code_type) then
          str[code_type] << [report.write_out()]
        else
          str[code_type] = []
          str[code_type] << [report.write_out()]
        end }
      
      # Write everything to disc
      str.each{ |key, value|
        myFile = File.open("#{filename}.yuverif-#{"%03d" % key}.data",'w')
        myFile.write(str_header.join("\n"))
        myFile.write(value.join("\n"))
        myFile.close
      }
      
    end
    
    def datatypes
      datatypes = Hash.new
      @reports.each { |report|
        newtype = report.repHeader.code_type_str
        if datatypes.has_key?(newtype) then
          datatypes[newtype] += 1
        else
          datatypes[newtype] = 1
        end
      }
      return datatypes
    end

    def summary
      datatypes = self.datatypes
      str  = "\n SUMMARY\n"
      str += "-"*38 + "\n"
      datatypes.each {|key,value|
    	str += " #{key.ljust(21)}: \t#{value}\n" }
      str += "-"*38 + "\n"
      str += " Total number of reports: \t#{self.length}\n\n"
      return str
    end
    
    def stations
      stations = Hash.new
      @reports.each { |report|
        newtype = report.repHeader.station_id
        if stations.has_key?(newtype) then
          stations[newtype] += 1
        else
          stations[newtype] = 1
        end
      }
      return stations
    end  


    def to_dfile(outfile)
      bulletins   = @reports.map{|r| r.to_bulletin}
      messages    = RubyBUFRTools::OBSMessages.new(bulletins)
      pp messages
      messages.to_dfile(outfile,"DUMMY")
    end

    
    private
    
    def parse(str)
      # make something with str
      
      # extract info from str
      strbeg = str.index(HEADER)
      if strbeg==nil 
        raise RuntimeError,"error: no leading \"VOF: [...]\" in #{$filename}\n" 
      end


      strend = str.index("YUVERIF ends")
      strend = -1 if strend.nil?
      
      # remove TRAILER lines (sometimes there is a TRAILER in the middle
      # of the text and we need to take this out....)
      a = str[strbeg..strend].split(/#{TRAILER}.*0\n/)
      a = a.delete_if{|x| x.split(/\n/).length == 1}
      #puts a.join("\n"+"="*50+"\n")
      str = a.join("\n")
      
      # now we can start parsing...
      number_of_runs = if str =~ Regexp.new(NUMBEROFRUNS + ".*$") then
                         $&.split(":")[1].to_i
                       else
                         1
                       end
      
      str = str.split("\n")
      
      ### extract the start date
      str[1] =~ /(\d+)/
      ts = $1 
      ts = ts.unpack("a4a2a2a2")
      basetime = Time.utc(ts[0],ts[1],ts[2],ts[3])
      verif_start,verif_end = str[2].scan(/.*\:(.*)\,.*\:(.*)/).flatten.collect{|i| i.to_f}
      
      # model run list
      strruntype = str.index(ModelRunTypes) + 1
      iline = strruntype + number_of_runs
      runs_info = str[strruntype..iline-1]
      
      iline = strruntype + number_of_runs
      txt_info = str.slice!(0..iline-1)
      
      runs_type  = []
      runs_start = []
      runs_endtime = []
      runs_mesh_width = []
      runs_nlevels  = []
      runs_description = []
      runs_info.each{ | currentline | 
    	items = currentline.split(" ")
        runs_type.push(items[1].to_i)
        strtime = items[2].unpack("a4a2a2a2")
        runs_start.push(Time.utc(strtime[0],strtime[1],strtime[2],strtime[3]))
        runs_endtime.push(items[3].to_f)
        runs_mesh_width.push(items[4].to_f)
        runs_nlevels.push(items[5].to_i)
        runs_description.push(items[6])
      }
      
      @info = {"text" => txt_info, "runs" => runs_info, "basetime" => basetime,
               "start" =>verif_start, "end" => verif_end}
      @runs = {"type" => runs_type, "start" => runs_start, "endtime" => runs_endtime,
               "mesh_width" => runs_mesh_width, "nlevels" => runs_nlevels,
               "description" => runs_description }
      
      str.shift # remove one empty line
      
      
      ### split all data lines in data items and remove empty lines
      #str.collect! { |theline|
      #  #splitted = theline.gsub("-"," -").split
      #  splitted = theline.split
      #  if splitted.length > 0 then
      #    splitted
      #  else
      #    nil
      #  end
      #}.compact!
      
      # remove empty lines
      str = str.select{|theline| theline =~ /\S/}

      # make reports from str
      # splits each report and creates Report objects
      iline = 0
      nlines = 1
      while not str.nil?

        if iline >= str.length-1
          break
        end
        theline = str[0].split
#         puts "theline"
#         p theline
        ### report header has 15 entries, otherwise reject
        if theline.length != 15 then
          next
        else
        end
        
        report_type = theline[1-1].to_i
        nlevels = [1,report_type].max
        report_type = [1,report_type].min	# multi-level report are type 1
        

        expected_columns = REPORT_TYPE_EXPECTED_COLUMNS[report_type][1..3] 
        
        nlines = 1 		# report header
        nlines += nlevels 	#  +report body
        nlines += nlevels * (number_of_runs.to_f/expected_columns[2]).ceil # optional rep. body
        # The following lines adds one line per block to the number of lines in 
        # the report if it is a multi-level report (the have on line with 9999 at
        # the end
        nlines += [0,report_type].max * (number_of_runs.to_f/expected_columns[2]).ceil
        
        nlines +=1 if report_type == 0 # SYNOP reports: reg report body is on two lines 
        
        # Now we should have the lines belonging to the next report
        strReport = str.slice!(0..nlines-1)

        if report_type == 0 then # SYNOP report store regular report body on two lines
          strReport[1] = strReport[1]+strReport[2]
          strReport.delete_at(2)
        end

        #puts nlevels
        #puts strReport.inspect
        #puts strReport[1..nlevels].inspect
        # Regular report body can have some formatting problem
        # The fortran based format looks like that:
        #	Line 1: '([I5,I5,I5,I5,I7,I6,I3,I2,I11,I4,I4,I2,I5,I11,I11,I5])'
        #	Line 2: '([I5,I5,I5,I4,I4,I5])' only for SYNOP reports
        format = [5,5,5,5,7,6,3,2,11,4,4,2,5,11,11,5,5,5,5,4,4,5]
        # here we process the first line(s) of the regular report body
        (1..nlevels).each{|i|
          thisline = strReport[i]
          #strReport[1..nlevels].collect!{|theline|
          #theline = theline.join(" ")
          #puts 
          regexp = Regexp.new(format[0..expected_columns[0]-1].collect{ |elem| "("+"."*elem+")"}.join)
          #puts regexp
          #puts theline.inspect
          #dummy = theline.scan(regexp).flatten.collect{|elem| elem.to_i}
          #puts dummy.inspect
          #puts dummy.class
          strReport[i] = thisline.scan(regexp).flatten.join(' ')
        }
        
        #if report_type == 0 then # SYNOP report store regular report body on two lines
        #  regexp = Regexp.new(format[1].collect{ |elem| "("+"."*elem+")"}.join)
        #  puts regexp
        #  #puts strReport[2].inspect
        #  #theline = strReport[2].join(" ")
        #  puts theline.inspect
        #  exit
        #  strReport[2] = strReport[2].join(" ").scan(regexp).flatten.collect{|elem| elem.to_i}
        #  
        #end
        
        strReport.collect! { |thisline| 
          thisline.gsub("-"," -").split
        }
        
        #puts "report.new with",strReport.inspect
        @reports.push(Report.new(strReport, self))
        
      end # while true

    end # method parse
    
  end # class YuVerif
  
  #
  # Reports found in the Yuverif file
  #
  class Report

    @@seed = 0

    attr_reader 	:repHeader, :basetime, :runs, :multi_level, :single_level, :report_type, \
    :regRepBody, :optRepBody, :nlevels,:missing
    #introduce hash of optional report bodies with key: initial time of run and value: Optional report body of this run
    
    
    def initialize(str, yuverif)

      #puts str.inspect
      #exit

      if str.class == String then
    	str = str.split("\n")
    	str[0] = str[0].split
      end
      
      report_type = str[0][0].to_i
      nlevels = [1,report_type].max
      report_type = [1,report_type].min	# multi-level report are type 1
      expected_columns = REPORT_TYPE_EXPECTED_COLUMNS[report_type][1..3] 
      
      @missing = MISSING
      
      #report header
      strhead = str.shift
      #    $stdout.print "Adding report with header: #{strhead.inspect}; "
      #    $stdout.print(".") if yuverif.reports.length.modulo(10) == 0
      
      @basetime   = yuverif.basetime
      @runs = yuverif.runs["start"]
      run_endleadtimes = yuverif.runs["endtime"]
      @nlevels 	= nlevels
      @multi_level  = (report_type == 1)
      @single_level = (report_type != 1)
      @report_type = report_type
      number_of_runs = @runs.length
      
      #puts strhead.inspect
      @repHeader    = ReportHeader.new(strhead, self)
      
      #extract regular report body from str
      #puts nlevels
      #puts str.inspect
      strregrep = str.slice!(0..nlevels-1)

      #puts strregrep.inspect
      #exit
      @regRepBody = RegReportBody.new(strregrep, self)
      
      #extract optional report body from str
      @optRepBody = {}
      # fill @optRepBody
      
      
      dummy = []
      while not str.length == 0
        #puts "-------",str.inspect
        #puts "==========="
        #puts dummy.inspect
        str.slice!(0..nlevels-1).transpose.each{|elem| dummy << elem}
        str.shift if report_type == 1 # multi-level reports have an additional line
      end
      
      #puts "DUMMY:\n #{dummy.inspect}"# if report_type == 0 

      # initialize number_of_runs instances of OptReportBody objects.
      (0..number_of_runs-1).each{|i|
        therun = @runs[i]
        theendleadtime = run_endleadtimes[i]
        lead_time = theendleadtime*3600 - (yuverif.info['end'].to_i*3600-@repHeader.obs_time_rel*60)
        #puts "run start time: #{therun}"
        #puts "obs_time_rel  : #{@repHeader.obs_time_rel}"
        #puts "lead_time     : #{lead_time/3600}"
        @optRepBody[therun.strftime("%Y%m%d%H")] = 
                        OptReportBody.new(dummy.slice!(0..expected_columns[1]-1),self,lead_time)
      }
    end
    
    def time
      @repHeader.time
    end
    
    def code_type
      @repHeader.code_type
    end
 
    def obs_status
      @repHeader.obs_status
    end

    def code_type_str
      @repHeader.code_type_str
    end
    
    def write_txt(format='old')
      str = [@repHeader.write_txt(format)]
      str << @regRepBody.write_txt(self)
      @optRepBody.each{|key,val|
        str << val.write_txt(self) }
      
      str = str.join("\n")
      return str 
    end

    def write_csv_header()
      str = [@repHeader.write_csv_header()," lead_time"] 
      #str.concat([@regRepBody.write_csv_header(self)])
      #key = @optRepBody.keys[0]
      #val = @optRepBody[key]
      #str.concat([val.write_csv_header(self)])
      return str.join(',')
    end


    def write_csv(run=nil)
      if run.nil?
        str = [@repHeader.write_csv()]
      else
        str = [@repHeader.write_csv(),@regRepBody.write_csv(self),@optRepBody['%c'.stimef(run.utc)].write_csv(self)]
      end
      return str.join(',')
    end


    def write_out()
      format='new'
      str = [@repHeader.write_out()]
      str << @regRepBody.write_txt(self)
      @optRepBody.each{|key,val|
        str << val.write_txt(self) }
      
      str = str.join("\n")
      return str 
    end
    
    def print()
      puts "============================"
      puts "Report Header"
      puts "============================"
      @repHeader.print
      puts "----------------------------"
      puts "Report Body"
      puts "----------------------------"
      @regRepBody.print
      if @optRepBody.length > 0 then
        puts "----------------------------"
        puts "Increments"
        puts "----------------------------"
        @optRepBody.each{ |key,vals|
          puts "run: #{key}"
          vals.print }
      end
    end
  
    def to_bulletin()
      bulletin = {}
      if @repHeader.code_type_str =~ /AIRCRAFT/
        aircraft_chars = (48..57).to_a.concat((65..90).to_a).concat((97..122).to_a)
        n_chars = aircraft_chars.length
        r = Random.new(@@seed)
        s=''
        6.times{
          s+= aircraft_chars[r.rand(0..n_chars-1)].chr
        }
        bulletin["refsite"]         = s
        @@seed += 1
      else
        bulletin["refsite"]         = @repHeader.station_id
      end
      bulletin["reftime"]         = @repHeader.time
      bulletin["stationheight"]   = @repHeader.sta_alt
      bulletin["position"]        = Met::LatLon.new(@repHeader.sta_lat,@repHeader.sta_lon)
      bulletin["origcenter"]      = 78
      bulletin["nlevels"]         = 1
      bulletin["unitDUMMY"]       = "dum"
      bulletin["nameDUMMY"]       = "DUMMY"
      bulletin["DUMMY"]           = 0.0
      bulletin["missing_value"]   = -9999.0
      return bulletin
    end

  
  end # class Report
  
  class ReportHeader
    
    attr_reader :report_type,:station_id,:sta_lon,:sta_lat,:obs_time_rel,\
    :sta_alt,:mod_alt,:obs_type,:code_type,:code_type_str,:stat_char,:rep_flag,\
    :obs_status,:thres_qc,:mod_i,:mod_j,:time
    
    def initialize(str, report)
      if str.class == String then
    	item = str.split(" ")
      else
    	item = str
      end
      
      @report_type   = item[ 1-1].to_i
      @station_id    = item[ 2-1].to_s 		# sometimes there is text
      @sta_lon       = item[ 3-1].to_f/100	# [°]
      @sta_lat       = item[ 4-1].to_f/100	# [°]
      @obs_time_rel  = item[ 5-1].to_i		# [min]
      @sta_alt       = item[ 6-1].to_i		# [m]
      @mod_alt       = item[ 7-1].to_i		# [m]
      @obs_type      = item[ 8-1].to_i
      @code_type     = item[ 9-1].to_i
      @code_type_str = YU::CTYPS[@code_type]
      @stat_char     = item[10-1].to_i
      @rep_flag      = item[11-1].to_i
      @obs_status    = item[12-1].to_i
      @thres_qc      = item[13-1].to_i
      @mod_i         = item[14-1].to_i
      @mod_j         = item[15-1].to_i
      
      @time		   = report.basetime + @obs_time_rel * 60
      
    end
    
    def write_txt(format='old')
      
      #    PRINTORDER = {
      #    0 => "id",
      #    1 => "lat",
      #    2 => "lon",
      #    3 => "reltime",
      #    4 => "time",
      #    5 => "sta_height",
      #    6 => "mod_height",
      #    7 => "obs_type",
      #    8 => "cod_type",
      #    9 => "type",
      #    10 => "status",
      #    11 => "report_flags",
      #    12 => "station_char",
      #    13=> "data"
      #  } 
      
      str = []
      str.concat([ @sta_lat,
                  @sta_lon,
                  @obs_time_rel ])
      str.concat([ if format != 'old' then	    
                     @time.utc.strftime("%Y%m%d%H%M")
                   else
                     @time.to_s
                   end])
      str.concat(
                 [@sta_alt,
                  @mod_alt,
                  @obs_type,
                  @code_type,
                  @report_type] )
      str.concat( [
                    @obs_status,
                    ("%030b" % @stat_char.to_i).reverse.scan(/(......)/).flatten.join(":"),
                    ("%030b" % @rep_flag.to_i).reverse.scan(/(......)/).flatten.join(":"),
                    @thres_qc,
                    @mod_i, @mod_j
                  ] ) if format != 'old'
      
      str = [@station_id.to_s].concat([str.join(",")])
    end


    def write_csv_header()
      str =      ["        time", "station_id", "sta_lon", "sta_lat", "sta_alt"]
      str.concat(%w[mod_i mod_j mod_alt obs_type code_type report_type obs_status])
      str.concat(["                         stat_char", "                          rep_flag", "thres_qc"])
      return str.join(",")
    end

    def write_csv()

      str = []
      str.concat([@time.utc.strftime("%Y%m%d%H%M"),
                  "%10d"    % @station_id.to_i,
                  "%7.2f"   % @sta_lon.to_f,
                  "%7.2f"   % @sta_lat.to_f])
      str.concat(
                 ["%7d"     % @sta_alt.to_i,
                  "%5d"     % @mod_i.to_i, 
                  "%5d"     % @mod_j.to_i,
                  "%7d"     % @mod_alt.to_i,
                  "%8d"     % @obs_type.to_i,
                  "%9d"     % @code_type.to_i,
                  "%11d"    % @report_type.to_i] )
      str.concat( ["%10d"    % @obs_status.to_i,
                   ("%030b" % @stat_char.to_i).reverse.scan(/(......)/).flatten.join(":"),
                   ("%030b" % @rep_flag.to_i).reverse.scan(/(......)/).flatten.join(":"),
                   "%8d" % @thres_qc.to_i])
      return str.join(",")
    end


    def write_out()
      
      #    PRINTORDER = {
      #    0 => "id",
      #    1 => "lat",
      #    2 => "lon",
      #    3 => "reltime",
      #    4 => "time",
      #    5 => "sta_height",
      #    6 => "mod_height",
      #    7 => "obs_type",
      #    8 => "cod_type",
      #    9 => "type",
      #    10 => "status",
      #    11 => "report_flags",
      #    12 => "station_char",
      #    13=> "data"
      #  } 
      
      str = []
      str.concat([@station_id,
                  @sta_lat,
                  @sta_lon,
                  @obs_time_rel ])
      str.concat([@time.utc.strftime("%Y%m%d%H%M")])
      str.concat(
                 [@sta_alt,
                  @mod_alt,
                  @obs_type,
                  @code_type,
                  @report_type] )
      str.concat( [
                    @obs_status,
                    ("%030b" % @stat_char.to_i).reverse.scan(/(......)/).flatten.join(":"),
                    ("%030b" % @rep_flag.to_i).reverse.scan(/(......)/).flatten.join(":"),
                    @thres_qc,
                    @mod_i, @mod_j
                  ] )
      
      str = str.join(",")
    end

    
    def print()
      puts self.inspect
    end
    
  end
  
  class RegReportBody
    
    attr_reader :u, :v, :temp, :rh, :pressure , :height, :finite_obs_err_flag, \
    :thres_qc_flag, :main_flag, :level_id, :pressure_code, :pressure_tendency, \
    :low_cloud_cover, :horiz_visibility, :combined_cloud, :cwpt_flag, \
    :precipitation_amount, :min_temp, :max_temp, :min_ground_temp, :max_wind_gust, \
    :max_min_speed, :global_radiation, :active_obs
    
    
    def initialize(str, report)
      
      #puts str.class
      #puts str.inspect

      report_type = report.report_type
      @report = report
      @u     = []
      @v     = []
      @temp  = []
      @rh    = []
      @pressure         = []
      @height           = []
      @finite_obs_err_flag = []
      @thres_qc_flag    = []
      @main_flag        = []
      @level_id         = []
      @pressure_code    = []
      @pressure_tendency = []
      @low_cloud_cover  = []
      @horiz_visibility = []
      @combined_cloud   = []
      # flag word for cloud, weather, precipitation, extreme temp
      @cwpt_flag        = [] 
      @precipitation_amount = []
      @min_temp         = []
      @max_temp         = []
      @min_ground_temp  = []
      @max_wind_gust    = []
      @max_min_speed    = []
      @global_radiation = []
      # GPS variables:
      @derived_IWV      = []
      @reported_IWV     = []
      @zenith_wet_delay = []
      @active_obs       = []
      if str.class == String then
    	str = str.split("\n").collect{|elem|
          elem.gsub!("-"," -").split }
	# one should add .transpose and use the same technique
	# as for OptRepBody to push data into the variables
      end
      
      str.each{ |line|
	line = line.split if str.class == String
	
	if report_type == 4 then # GPS report are special
	  @derived_IWV.push(Yuverif.convert(line[1-1],0.01,-9))
          @reported_IWV.push(Yuverif.convert(line[2-1],0.01,-9))
        else
          #puts line.class
          #puts line.inspect
          #puts line[0]
          #puts line[1]
          #exit
          @u.push(Yuverif.convert(line[1-1],0.1))
          @v.push(Yuverif.convert(line[2-1],0.1))
        end
        @temp.push(Yuverif.convert(line[3-1],0.1,-9))
        @rh.push(Yuverif.convert(line[4-1],0.1,-9))
        @pressure.push(Yuverif.convert(line[5-1],1,-9))
        @height.push(Yuverif.convert(line[6-1],1,-999))
        @finite_obs_err_flag.push(Yuverif.convert(line[7-1],1,-9))
        @thres_qc_flag.push(Yuverif.convert(line[8-1],1,-9))
        @main_flag.push(Yuverif.convert(line[9-1],1,-9))
        if report_type == 0 then
          @pressure_code.push(Yuverif.convert(line[10-1],1,-9))
        else
          @level_id.push(Yuverif.convert(line[10-1],1,-9))
        end
        if report_type < 1 then # the following lines do no apply to multi-level reports
          if  report_type == 4 then # GPS report are special
            @zenith_wet_delay.push(Yuverif.convert(line[11-1],1,-9))
          else
            @pressure_tendency.push(Yuverif.convert((line[11-1].to_i-500).to_s,0.1,-9)) # [Pa/3H]
          end
          if [0,-1].include?(report_type) then # short surface-level or SYNOP
            @low_cloud_cover.push(Yuverif.convert(line[12-1],1,-9))
            @horiz_visibility.push(Yuverif.convert(line[13-1],1,-9))
            @combined_cloud.push(Yuverif.convert(line[14-1],1,-9))
            @cwpt_flag.push(Yuverif.convert(line[15-1])) # flag word for cloud, weather, precipitation, extreme temp
            @precipitation_amount.push(Yuverif.convert(line[16-1],0.1,-9))
            if report_type == 0 then # SYNOP only
              @min_temp.push(Yuverif.convert(line[17-1],0.1,-9))
              @max_temp.push(Yuverif.convert(line[18-1],0.1,-9))
              @min_ground_temp.push(Yuverif.convert(line[19-1],0.1,-9))
              @max_wind_gust.push(Yuverif.convert(line[20-1]))
              @max_min_speed.push(Yuverif.convert(line[21-1]))
              @global_radiation.push(Yuverif.convert(line[22-1],10000,-9)) # [J/m^2]
            end
          end
        end
      }
      
    end
    
    def no_active_obs?
      determine_active_obs if @active_obs.empty?
      a = @active_obs.select{|x| x.join =~ /[A-Z]+/}.empty?
      #puts "no active obs (equivalent to obs_status = 2)" if a
      return a
    end


    def determine_active_obs

      return unless @active_obs.empty?
      # loop over all levels
      @finite_obs_err_flag.each_index{|i|
        wind,temp,hum,pg = "-","-","-","-"
        # if report is active...
        if @report.obs_status == 0 
          fin_obs_err_flag = ('%06b'%@finite_obs_err_flag[i]).reverse.scan(/(.)(.)(.)(.)(.)(.)/).flatten
          thres_qc_flag    = ('%04b'%@thres_qc_flag[i]).reverse.scan(/(.)(.)(.)(.)/).flatten
          wind = "W" if (fin_obs_err_flag[0].to_i == 1 and thres_qc_flag[0].to_i == 0)
          temp = "T" if (fin_obs_err_flag[1].to_i == 1 and thres_qc_flag[1].to_i == 0)  
          hum  = "H" if (fin_obs_err_flag[2].to_i == 1 and thres_qc_flag[2].to_i == 0)
          if (@report.report_type <= 0 and @report.report_type >= -2)
            pg = "P" if (fin_obs_err_flag[3].to_i == 1 and thres_qc_flag[3].to_i == 0)
          else
            pg = "G" if (fin_obs_err_flag[3].to_i == 1 and thres_qc_flag[3].to_i == 0)
          end
        end
        @active_obs.push([wind,temp,hum,pg])
      }

      

    end



    def write_csv_header(report)
      report_type = report.report_type
      str =      %w[u v temp rh pressure height finite_obs_err_flag thres_qc_flag main_flag level_id]
      str.concat(%w[pressure_code pressure_tend low_cloud_cover horiz_vis comb_cloud cwpt_flag prec_amount])
      str.concat(%w[min_temp max_temp min_ground_temp max_wind_gust max_min_speed global_rad])
      if report_type == -4 
	# GPS variables:
    	str[1-1] = "derived_IWV"
    	str[2-1] = "reported_IWV"
    	str[11-1] = "zenith_wet_delay"
      end
      return str.join(",")
    end

    def write_csv(report)
      report_type = report.report_type
      
      str  = [@u,
              @v,
              @temp,
              @rh,
              @pressure,
              @height,
              @finite_obs_err_flag,
              @thres_qc_flag,
              @main_flag,
              @level_id,
              @pressure_code,
              @pressure_tendency,
              @low_cloud_cover,
              @horiz_visibility,
              @combined_cloud,
              @cwpt_flag, # flag word for cloud, weather, precipitation, extreme temp
              @precipitation_amount,
              @min_temp,
              @max_temp,
              @min_ground_temp,
              @max_wind_gust,
              @max_min_speed,
              @global_radiation]
      
      if report_type == -4 
	# GPS variables:
    	str[1-1] = @derived_IWV
    	str[2-1] = @reported_IWV
    	str[11-1] = @zenith_wet_delay
      end
      
      return str.join(",")
    end
    

    def write_txt(report)
      str = []
      
      report_type = report.report_type
      
      temp = [@u,
              @v,
              @temp,
              @rh,
              @pressure,
              @height,
              @finite_obs_err_flag,
              @thres_qc_flag,
              @main_flag,
              @level_id,
              @pressure_code,
              @pressure_tendency,
              @low_cloud_cover,
              @horiz_visibility,
              @combined_cloud,
              @cwpt_flag, # flag word for cloud, weather, precipitation, extreme temp
              @precipitation_amount,
              @min_temp,
              @max_temp,
              @min_ground_temp,
              @max_wind_gust,
              @max_min_speed,
              @global_radiation]
      
      if report_type == -4 
	# GPS variables:
    	temp[1-1] = @derived_IWV
    	temp[2-1] = @reported_IWV
    	temp[11-1] = @zenith_wet_delay
      end
      
      
      temp.collect!{|elem| if elem.length == 0 then nil else elem end}.compact!
      temp = temp.transpose
      str << temp.join(",")
      return str
    end

    def print()
      puts self.inspect
    end
    
  end # Class RegRepBody
  
  
  class OptReportBody
    
    attr_reader :iwv,:u,:v,:temp,:rh,:geopotential,:pressure,:lead_time,
        :total_cloud_cover,:low_cloud_cover,:horiz_visibility,:precipitation_amount,
        :min_temp,:max_temp,:min_ground_temp,:max_wind_speed_gust,:global_radiation 
    
    def initialize(str, report, lead_time)
      
      report_type = report.report_type
      @lead_time = lead_time
      #  @u = []
      #  @v = []
      #  @iwv = []
      #  @temp = []
      #  @rh = []
      #  @pressure = []
      #  @geopotential = []
      #  @total_cloud_cover = []
      #  @low_cloud_cover = []
      #  @horiz_visibility = []
      #  @precipitation_amount = []
      #  @min_temp = []
      #  @max_temp = []
      #  @min_ground_temp = []
      #  @max_wind_speed_gust = []
      #  @global_radiation = []
      
      if str.class == String then
    	str = str.split("\n").collect{|elem|
          elem.gsub!("-"," -").split }.transpose
      end
      
      if report_type == 4 then # GPS report are special
        @iwv = Yuverif.strarrtonumarr(str[1-1],0.01,-9)
      else
        @u = Yuverif.strarrtonumarr(str[1-1],0.01,9999)
        @v = Yuverif.strarrtonumarr(str[2-1],0.01,9999)
        @temp = Yuverif.strarrtonumarr(str[3-1],0.01,9999)
        @rh = Yuverif.strarrtonumarr(str[4-1],0.1,9999)
        if report_type != -3 then
          if report_type == 1 then
            @geopotential = Yuverif.strarrtonumarr(str[5-1],1,9999)
          else
            @pressure = Yuverif.strarrtonumarr(str[5-1],1,9999)
          end
          if report_type < 1 and report_type != -2 then
            @total_cloud_cover = Yuverif.strarrtonumarr(str[6-1],1.0,-9)
            @low_cloud_cover = Yuverif.strarrtonumarr(str[7-1],1.0,-9)
            @horiz_visibility = Yuverif.strarrtonumarr(str[8-1],10,-9)
            @precipitation_amount = Yuverif.strarrtonumarr(str[9-1],0.1,-9)
            if report_type == 0 then
              @min_temp = Yuverif.strarrtonumarr(str[10-1],0.1,9999)
              @max_temp = Yuverif.strarrtonumarr(str[11-1],0.1,9999)
              @min_ground_temp = Yuverif.strarrtonumarr(str[12-1],0.1,9999)
              @max_wind_speed_gust = Yuverif.strarrtonumarr(str[13-1],1.0,9999)
              @global_radiation = Yuverif.strarrtonumarr(str[14-1],10000,9999)
            end
          end
        end
      end
    end
    
    def write_txt(report)
      str = []
      
      report_type = report.report_type
      
      str.concat( [@u,
                   @v,
                   @temp,
                   @rh,
                   @pressure,
                   @total_cloud_cover,
                   @low_cloud_cover,
                   @horiz_visibility,
                   @precipitation_amount,
                   @min_temp,
                   @max_temp,
                   @min_ground_temp,
                   @max_wind_speed_gust,
                   @global_radiation] )
      
      # Upper-air, multi-level
      str[5-1] = @geopotential if report_type == 1
      # GPS
      str[1-1] = @iwv if report_type == -4
      
      str.compact!
      str = str.transpose.join(",")
      return str
    end


    def write_csv_header(report)
      report_type = report.report_type
      str = %w[u v temp rh pressure total_cloud_cover low_cloud_cover]
      str.concat(%w[horiz_visibility precipitation_amount min_temp max_temp min_ground_temp max_wind_speed_gust global_radiation])
      # Upper-air, multi-level
      str[5-1] = "geopotential" if report_type == 1
      # GPS
      str[1-1] = "iwv" if report_type == -4 

      str.concat(["lead_time"])

      return str.join(",")
    end



    def write_csv(report)
      str = []
      
      report_type = report.report_type
      
      str.concat( [@u,
                   @v,
                   @temp,
                   @rh,
                   @pressure,
                   @total_cloud_cover,
                   @low_cloud_cover,
                   @horiz_visibility,
                   @precipitation_amount,
                   @min_temp,
                   @max_temp,
                   @min_ground_temp,
                   @max_wind_speed_gust,
                   @global_radiation] )
      
      # Upper-air, multi-level
      str[5-1] = @geopotential if report_type == 1
      # GPS
      str[1-1] = @iwv if report_type == -4
      
      
      str.concat([@lead_time])
      str = str.join(",")
      return str
    end    



    def print
      puts self.inspect
    end
    
  end # Class OptRepBody
  

  class BlackList

    attr_reader :raw_list, :blacklist, :whitelist

    def initialize(list)
      if (File.file?(list)) then
        @raw_list = File.open(list,"r").readlines
      elsif (list.is_a?(String)) then
        @raw_list = list.split("\n")
      elsif (list.is_a?(Array)) then
        @raw_list = list
      else
        raise ArgumentError,"Error: list needs to be a file, a string or an array"
      end
      
      parse
      
    end

    private

    def parse(list=@raw_list)
      
      ind_whitelist = list.rindex("WHITELIST\n")
      
      last_black = ind_whitelist.nil? ? list.length : ind_whitelist

      # parse blacklist
      @blacklist = []
      list[1..last_black-1].each{|line|
        a = line.split
        hash = {}
        hash["station_id"]   = a[0]
        hash["obs_type"]     = a[1]
        hash["limit_geopot"] = a[2..3]
        hash["limit_wind"]   = a[4..5]
        hash["limit_temp"]   = a[6..7]
        hash["limit_hum"]    = a[8..9]
        @blacklist.push(hash)
      }

      return if ind_whitelist.nil?
      
      # parse whitelist
      @whitelist = []
      list[last_black+1..-1].each{|line|
        a = line.split
        hash = {}
        hash["station_id"]   = a[0]
        hash["obs_type"]     = a[1]
        hash["obs_code_type"]= a[2]
        @whitelist.push(hash)
      }
        
    end
    

    
  end # class BlackList





end # module
