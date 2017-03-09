package utils

import java.nio.file.Paths

import com.typesafe.config.ConfigFactory

/**
  * Created by jdelafon on 05/07/16.
  */
object Constants {

  private def getConfigPath(key:String): String = {
    Paths.get(ConfigFactory.load().getString(key)).toAbsolutePath.toString
  }

  val TEST_MODE: Boolean = ConfigFactory.load().getBoolean("env.TEST_MODE")
  val BAM_PATH: String = getConfigPath("env.BAM_PATH")
  val TEMP_BAM_DIR: String = getConfigPath("env.TEMP_BAM_DIR")
  val BAM_BAI_REGEX = ".*?\\.bam.*"
  val KEY_REGEX = """^([\w\d:.-]+?)(.bam){0,1}(.bai|.idx){0,1}$""".r

}
