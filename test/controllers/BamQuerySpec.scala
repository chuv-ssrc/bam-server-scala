package controllers

import javax.inject.Inject

import controllers.generic.BamQueryController
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Check error handling by the generic BamQueryCntroller.
  * Use IndexController for that.
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

    // See IndexSpec
//    "fail if the key is not found in database (POST)" in {
//      val body: JsValue = Json.parse(s"""{"key": "xxx"}""")
//      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
//      status(response) mustBe INTERNAL_SERVER_ERROR
//    }
//
//    "fail if the corresponding BAM file is not found on disk (POST)" in {
//      val body: JsValue = Json.parse(s"""{"key": "notherekey"}""")
//      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
//      status(response) mustBe INTERNAL_SERVER_ERROR
//    }


  }

}
