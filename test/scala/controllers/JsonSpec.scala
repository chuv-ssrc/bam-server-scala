package controllers

import org.scalatest.BeforeAndAfter
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Test JsonController.
  * The test bam has coordinates in range chr1:761997-762551.
  */
class JsonSpec extends PlaySpec with OneAppPerSuite with BeforeAndAfter {

  val token = "asdf"
  val testkey = "testkey"
  val body: JsValue = Json.parse(s"""{"key": "$testkey"}""")
  val headers = (AUTHORIZATION -> s"Bearer $token")

  "JsonController" should {

    "provide reads in JSON format if a region is given (POST)" in {
      val request = FakeRequest(POST, "/bam/json?region=chr1:761997-762551").withJsonBody(body).withHeaders(headers)
      val response = route(app, request).get
      status(response) mustBe OK
      //contentAsBytes(response).length must be > 100000
    }

    "provide reads in JSON format if a region is given (GET)" in {
      val request = FakeRequest(GET, s"/bam/json/$testkey/$token?region=chr1:761997-762551")
      val response = route(app, request).get
      status(response) mustBe OK
      //contentAsBytes(response).length must be > 100000
    }

  }

}