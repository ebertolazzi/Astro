#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"

FileUtils.rm_rf   "bin"
FileUtils.mkdir_p "bin"

FileUtils.cp_r "../src/.", "./src";
FileUtils.cp_r "../submodules/GenericContainer/src/.", "./src/GenericContainer";
FileUtils.cp_r "../submodules/GenericContainer/include/.", "./src";
FileUtils.cp_r "../submodules/UtilsLite/src/.", "./src";
FileUtils.cp_r "../submodules/quarticRootsFlocke/src/.", "./src";

FileUtils.cp "../license.txt", "license.txt"
