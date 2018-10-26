var kafka = require('kafka-node'),
    Producer = kafka.Producer,
    KeyedMessage = kafka.KeyedMessage,
    client = new kafka.Client(),
    producer = new Producer(client),
    km = new KeyedMessage('key', 'message ' + new Date().toDateString()),
    payloads = [
        { topic: 'topic2', messages: JSON.stringify({type:'S', name: 'Galileo', date:new Date(), attrs:{runner:[{name:'h1', track:'galway'},{name:'h2', track:'sligo'}]}}), partition: 0 , attributes : 2}
        // ,{ topic: 'topic2', messages: ['hello ' + new Date().toDateString(), 'world ' + new Date().toDateString(), km] , attributes : 2}
        // ,{ topic: 'topic3', messages: [new KeyedMessage({type:'S', name: 'Galileo', date:new Date()}, {runners:[{name:'h1', track:'galway'},{name:'h2', track:'sligo'}]})] }
    ];

producer.on('ready', function () {
    producer.send(payloads, function (err, data) {
        console.log(data);
    });
});

producer.on('error', function (err) {})