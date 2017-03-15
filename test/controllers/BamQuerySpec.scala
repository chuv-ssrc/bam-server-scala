package controllers

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Check error handling by the generic BamQueryCntroller.
  */
class BamQuerySpec extends PlaySpec with OneAppPerSuite {

  "`keyToBamRequest`" should {

    "return a BamRequest for the given sample key" in {

    }

    "complain if the key is not found in database" in {

    }

    "complain if the corresponding BAM file is not found on disk" in {

    }

  }

  "`parseBamRequestFromPost`" should {

    "return a BamRequest from the JSON body of a request" in {

    }

    "complain if there is no JSON body" in {

    }

    "complain if there is no 'key' in the body" in {

    }

  }

}
