# SSHing into the repository

This is straightforward if you follow the following guide for generating an SSH keypair.
https://docs.gitlab.com/ee/ssh/#ssh-on-the-gitlab-server

ie. `ssh-keygen -t rsa -b 4096 -C "email@example.com"`, where the email is optional.

The private key will be placed in your local `.ssh` dir, while the contents of the `.pub` file should be pasted into the SSH key section in your GitLab user account (ie. via the GUI).

Test that everything is working correctly by running: `ssh -Tv git@git.chalmers.se` verbatim.

You should receive a message like 
>Welcome to GitLab, @croft!
if everything is working correctly.


## Troubleshooting
If you are prompted for a password, the keys are not set up correctly.

Try running `ssh -Tvvv git@git.chalmers.se` (more v's means more verbosity...) to trace the point of failure. 


If you have multiple keys in your `.ssh` then the ssh will try each of them in turn. If too many fail (four, I believe) the git server will reject any further keys and will ask for a password instead.

You can try to
1) Remove some other keys.
2) Specify the key somehow (in ssh it is `ssh -i <privkey> git@etc`, I don't know about `git`)
3) Set up an ssh config file to associate the url with a specific key. Recommended. See https://superuser.com/a/232406

##### Setting up `config`

Create a file named  `config` in `~/.ssh` and give it restrictive permissions (ie. 600).

Enter something like the following. Modify to suit.
```
Host git.chalmers.se
        HostName git.chalmers.se
        User git
        IdentityFile ~/.ssh/gitlab_chalmers
```

Finally, restart SSH.

