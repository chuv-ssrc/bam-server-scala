package scala

import org.scalatestplus.play._
import play.api.test._
import play.api.test.Helpers._


class BasicSpec extends PlaySpec with OneAppPerTest {

  "Routes" should {

    "send 404 on a bad request" in  {
      route(app, FakeRequest(GET, "/fake")).map(status) mustBe Some(NOT_FOUND)
    }

  }

  "HomeController" should {

    "render the index page" in {
      val response = route(app, FakeRequest(GET, "/")).get
      status(response) mustBe OK
      contentType(response) mustBe Some("text/plain")
      contentAsString(response) must include ("BAM server operational.")
    }

  }

}
