package controllers

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


class IndexSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"

  "IndexController" should {

    "provide the BAM index (POST) if everything is right" in {
      val header = ("Authorization", s"Bearer $token")
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe OK
    }

    "provide the BAM index (GET) if everything is right" in {
      val response = route(app, FakeRequest(GET, s"/bai/testkey/$token")).get
      status(response) mustBe OK
    }

    "fail if the key is not known to the database" in {
      val header = ("Authorization", s"Bearer $token")
      val body: JsValue = Json.parse(s"""{"key": "xxx"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
    }

    "fail if the key exists but the bam file is not found" in {
      val header = ("Authorization", s"Bearer $token")
      val body: JsValue = Json.parse(s"""{"key": "notherekey"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
    }

  }

}
