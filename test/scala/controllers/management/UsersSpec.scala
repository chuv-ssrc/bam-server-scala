package scala.controllers.management

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import utils.DbUtils._

import scala.setup.{WithDatabase, WithToken}


/**
  * Test UsersController.
  */
class UsersSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {

  val db = dbContext.db

  "Add a user" should {

    "add rows to the users table" in {
      val body = Json.parse("""
        { "users": [{"username": "A"}, {"username": "B"}] }
      """)
      val count0 = countRows(db, "users")
      val response = route(app, FakeAdminRequest(PUT, "/users").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users")
      count must equal(count0 + 2)
    }

    "fail and insert nothing if one of the users already exists" in {
      val body = Json.parse("""
        { "users": [{"username": "A"}, {"username": "B"}] }
      """)
      val count0 = countRows(db, "users")
      val response = route(app, FakeAdminRequest(PUT, "/users").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "users")
      count must equal(count0)
    }

  }

  "Delete a user" should {

    "remove rows from the users table" in {
      val body = Json.parse("""
        { "users": ["A","B"] }
      """)
      val count0 = countRows(db, "users")
      val response = route(app, FakeAdminRequest(DELETE, "/users").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users")
      count must equal(count0 - 2)
    }

    "fail and delete nothing if one of the users does not exist" in {
      val body = Json.parse("""
        { "users": ["A","C"] }
      """)
      val count0 = countRows(db, "users")
      val response = route(app, FakeAdminRequest(DELETE, "/users").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      val count = countRows(db, "users")
      count must equal(count0)
    }

  }

}
