package utils

import java.nio.file.Paths

import com.typesafe.config.ConfigFactory

/**
  * Created by jdelafon on 05/07/16.
  */
object Constants {

  private val config = ConfigFactory.load()

  private def getConfigPath(key: String): String = {
    Paths.get(config.getString(key)).toAbsolutePath.toString
  }

  val TEST_MODE: Boolean = config.getBoolean("env.TEST_MODE")
  val BAM_PATH: String = getConfigPath("env.BAM_PATH")
  val TOKEN_USER_KEY: String = config.getString("env.TOKEN_USER_KEY")
  val TOKEN_ISS_KEY: String = config.getString("env.TOKEN_ISS_KEY")

}
