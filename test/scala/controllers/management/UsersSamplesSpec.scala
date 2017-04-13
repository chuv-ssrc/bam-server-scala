package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test IndexController.
  */
class UsersSamplesSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  "Add an attribution" should {

    "add rows to the users_samples table" in {
      val body = Json.parse("""
        { "users_samples": [{"sample": "testkey", "username": "testuser"}, {"sample": "inactivekey", "username": "testuser"}] }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(PUT, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0 + 2)
    }

  }

  "Delete an attribution" should {

    "remove rows from the users_samples table" in {
      val body = Json.parse("""
        { "users_samples": [{"sample": "testkey", "username": "testuser"}, {"sample": "inactivekey", "username": "testuser"}] }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "users_samples")
      val response = route(app, FakeAdminRequest(DELETE, "/users_samples").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users_samples")
      count must equal(count0 - 2)
    }

  }

}
