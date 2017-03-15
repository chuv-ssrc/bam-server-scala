package controllers

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Test SamtoolsController.
  * The test bam has coordinates in range chr1:761997-762551.
  */
class SamtoolsSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"
  val headers = (AUTHORIZATION -> s"Bearer $token")

  "SamtoolsController" should {

    "Return everything if no region is given (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val headers = (AUTHORIZATION -> s"Bearer $token")
      val request = FakeRequest(POST, "/bam/samtools").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
    }

    "provide a slice of the BAM if a region is given (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val request = FakeRequest(POST, "/bam/samtools?region=chr1:0-10000").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
    }

//    "provide a slice of the BAM if everything is right (GET)" in {
//      val request = FakeRequest(GET, s"/bam/samtools/testkey/$token?range=0-65639")
//      val response = route(app, request).get
//      status(response) mustBe PARTIAL_CONTENT
//    }

  }

}