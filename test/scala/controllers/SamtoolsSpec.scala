package scala.controllers

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._
import utils.BamUtils.samtoolsExists


/**
  * Test SamtoolsController.
  * The test bam has coordinates in range chr1:761997-762551.
  */
class SamtoolsSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"
  val testkey = "testkey"
  val body: JsValue = Json.parse(s"""{"key": "$testkey"}""")
  val headers = (AUTHORIZATION -> s"Bearer $token")

  "SamtoolsController" should {

    "Return everything if no region is given (POST)" in {
      assume(samtoolsExists())
      val headers = (AUTHORIZATION -> s"Bearer $token")
      val request = FakeRequest(POST, "/bam/samtools").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
    }

    "provide a slice of the BAM if a region is given (POST)" in {
      assume(samtoolsExists())
      val request = FakeRequest(POST, "/bam/samtools?region=chr1:761997-762551").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
      contentAsBytes(response).length must be > 100000
    }

    "provide a slice of the BAM if a region is given (GET)" in {
      assume(samtoolsExists())
      val request = FakeRequest(GET, s"/bam/samtools/$testkey/$token?region=chr1:761997-762551")
      val response = route(app, request).get
      status(response) mustBe OK
      contentAsBytes(response).length must be > 100000
    }

    "return 0 bytes if the region is out of range (POST)" in {
      assume(samtoolsExists())
      val request = FakeRequest(POST, "/bam/samtools?region=chr1:0-10000").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
      contentAsBytes(response).length must equal(4645)
    }

    "return 0 bytes if the region is out of range (POST)(2)" in {
      assume(samtoolsExists())
      val request = FakeRequest(POST, "/bam/samtools?region=chr1:100000000-200000000").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
      contentAsBytes(response).length must equal(4645)
    }

  }

}