PARSE_APP_HOST=localhost
PARSE_APP_PORT=1337
CLOUD_DIR=/home/ian/Documents/sandbox/horsetracker-parse
APPLICATION_ID=test
MASTER_KEY=test
export NODE_PATH=/usr/local/lib/node_modules
sudo npm install -g parse-server mongodb-runner parse-dashboard parse
mongodb-runner start &
sleep 5
parse-server --appId $APPLICATION_ID --masterKey $MASTER_KEY --cloud $CLOUD_DIR/main.js --port $PARSE_APP_PORT & 
sleep 5
parse-dashboard --appId $APPLICATION_ID --masterKey $MASTER_KEY --serverURL "http://$PARSE_APP_HOST:$PARSE_APP_PORT/parse" --appName $APPLICATION_ID &
sleep 5


Parse.Cloud.define("hello", function(request, response) {
  response.success("Hello world, Ian!");
});

parse-dashboard --appId $APPLICATION_ID --masterKey $MASTER_KEY --serverURL "http://localhost:23740/parse" --appName demo &




curl -X POST \
  -H "X-Parse-Application-Id: test" \
  -H "Content-Type: application/json" \
  -d '{}' \
  http://$PARSE_APP_HOST:$PARSE_APP_PORT/parse/functions/hello
