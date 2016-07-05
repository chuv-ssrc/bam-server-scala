package utils

import java.nio.file.Paths

import com.typesafe.config.ConfigFactory

/**
  * Created by jdelafon on 05/07/16.
  */
object Constants {

  private def getConfigValue(key:String): String = {
    Paths.get(ConfigFactory.load().getString(key)).toAbsolutePath.toString
  }

  val BAM_PATH:String = getConfigValue("env.BAM_PATH")
  val APACHE_TEMP_BAM_DIR:String = getConfigValue("env.APACHE_TEMP_BAM_DIR")
  val APACHE_TEMP_BAM_URL:String = getConfigValue("env.APACHE_TEMP_BAM_URL")
  val TEMP_BAM_DIR:String = getConfigValue("env.TEMP_BAM_DIR")
  val BAM_BAI_REGEX = ".*?\\.bam.*"

}
