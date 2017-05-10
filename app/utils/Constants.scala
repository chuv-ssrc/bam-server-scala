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

}
