package controllers

import akka.stream.Materializer
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


class IndexSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"

  "IndexController" should {

    "provide the BAM index (POST) if everything is right" in {
      val header = ("Authorization", s"Bearer $token")
      val body: JsValue = Json.parse(s"""{"key": "aaaa"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe OK
    }

    "fail if the key is not known to the database" in {
      val header = ("Authorization", s"Bearer $token")
      val body: JsValue = Json.parse(s"""{"key": "aaaa"}""")
      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
    }


    "provide the BAM index (GET)" in {
      val response = route(app, FakeRequest(GET, "/bai/aaaa")).get
      status(response) mustBe OK
    }

  }

//  "HomeController" should {
//
//    "render the index page" in {
//      val response = route(app, FakeRequest(GET, "/")).get
//      status(response) mustBe OK
//      contentType(response) mustBe Some("text/plain")
//      contentAsString(response) must include ("BAM server operational.")
//    }
//
//  }


}
