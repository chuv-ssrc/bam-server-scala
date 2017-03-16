package scala.controllers

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Test IndexController.
  */
class IndexSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"
  val testkey = "testkey"
  val header = (AUTHORIZATION, s"Bearer $token")


  "`parseBamRequestFromPost`" should {

    "fail if there is no JSON body" in {
      val response = route(app, FakeRequest(POST, "/bai").withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No key found in request body")
    }

    "fail if there is no 'key' in the body" in {
      val body: JsValue = Json.parse(s"""{"yyyy": "xxx"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No key found in request body")
    }

  }


  "`keyToBamRequest`" should {

    "fail if the key is not known to the database (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "xxx"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found on disk (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "notherekey"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

    "fail if the key is not known to the database (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bai/xxx/$token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bai/notherekey/$token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

  }


  "IndexController" should {

    "provide the BAM index if everything is right (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "$testkey"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY  // application/octet-stream
    }

    "provide the BAM index if everything is right (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bai/$testkey/$token")).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

  }

}
