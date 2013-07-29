function hz=mel2hz(mel)
    hz=(exp(mel/1127.01048)-1)*700;
return