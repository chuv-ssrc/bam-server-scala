package scala

import org.scalatestplus.play._
import play.api.test._
import play.api.test.Helpers._
import java.nio.file.Paths


import org.scalatestplus.play._

import utils.Common._


class UtilsSpec extends PlaySpec {

  "`isOnDisk`" should {

    "Return whether the file is on disk" in {
      isOnDisk(Paths.get("asdf").toFile) must be (false)
      isOnDisk(Paths.get("scripts/onDisk.sh").toFile) must be (true)
    }

  }

}