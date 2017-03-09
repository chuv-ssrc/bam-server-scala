import org.scalatestplus.play._
import play.api.test._
import play.api.test.Helpers._
import play.api.db.Database
import java.nio.file.Paths


//import controllers.BamController

import org.scalatest.mock.MockitoSugar
import org.scalatestplus.play._


/**
  * Add your spec here.
  * You can mock out a whole application including requests, plugins etc.
  * For more information, consult the wiki.
  */
class UnitTestsSpec extends PlaySpec with MockitoSugar {

  "BamController" should {

    "Return whether the file is on disk" in {
//      val ctrl = new BamController(mock[Database])
//      ctrl.isOnDisk(Paths.get("asdf").toFile) must be (false)
//      ctrl.isOnDisk(Paths.get("scripts/onDisk.sh").toFile) must be (true)
    }

  }

}