package utils

import java.sql._
import com.googlecode.protobuf.format.JsonFormat
import protocol.Command.ResponseCommand
import scala.collection.JavaConversions._

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 29.10.13
 * Time: 9:19
 * To change this template use File | Settings | File Templates.
 */

object DbUtils {

  def saveWorkTask(work: Any, prep: PreparedStatement, conn: Connection): Int = {
    prep.setString(1, work.toString)
    prep.setString(2, "")
    prep.setString(3, "in progress")
    prep.executeUpdate

    val statement = conn.createStatement()
    val resultSet = statement.executeQuery("SELECT LAST_INSERT_ID()")

    if (resultSet.next()) {
      resultSet.getInt(1)
    } else -1
  }

  def updateTask(id: Int, result: Any, prep: PreparedStatement, conn: Connection) = {
    val builder = ResponseCommand.newBuilder()
    JsonFormat.merge(result.toString, builder)
    val responseCommand = builder.build()

    prep.setString(1, responseCommand.getDataList.map(JsonFormat.printToString).mkString(","))
    prep.setString(2, responseCommand.getStatus)
    prep.setString(3, id.toString)
    prep.executeUpdate
  }

  def findTask(id: Int, prep: PreparedStatement, conn: Connection) : String = {
    prep.setString(1, id.toString)
    val rs = prep.executeQuery()
    if (rs.next()) {
      rs.getString("status") match {
        case m @ "in progress" => m
        case m @ _ => "{\"status\": \"" + rs.getString("status") + "\", data: [" + rs.getString("result") + "]}"
      }
    } else {
      "Your job id not found"
    }
  }
}
