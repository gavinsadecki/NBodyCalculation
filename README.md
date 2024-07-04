# NBodyCalculation
This was a final project for one of my classes. I decided to revisit and clean up the code to show improvement of my skills over time.
The project was for numerical methods which discusses how computers handle numbers and using programming to compute answers to problems.
It was the final project for the class and it had asked us to compute the positions of planetary bodies using power series to estimate the terms.
I believe you could do it easier with kinematic equations but I haven't tried or looked further into it. I believe in class we went over this problem given 2 different
planetary bodies. The final project was to either expand the code to use 3 bodies or to make it support n bodies. The calculations I did by hand and saw a pattern to how you can generate the terms
given n different bodies and went for that route. I don't believe I have the notes anymore so I don't have the work I did to find the pattern but the code should be close enough.

## 2018
This was written when my experience was very informal and still in school. I knew how to use the language to get projects done but it did not look the best and might not be the most efficient. It was
very scrappy in how I have it graphing the output of the data. I left variables in there that I could easily swap in and out to get the output I need but this isn't a good way to do it.
For one, I don't remember how I had things organized or what I was changing in the past. There's also almost no documentation as to what any of this does or what the data is so it'd be hard
for someone to reuse this.

## 2022
After gaining work experience and with other programming languages, I decided to revisit this project to see how much my skills have changed since I originally worked on this.
One big improvement is I was able to simplify the calculations more. I saw a way it could be reduced further and implemented that here. It now looks way more readable.
While in school, I didn't use objects that often and never really saw the potential of it. I had lessons on it in Python but never did much beyond that. While working, I used C#
and I've grown to appreciate objects for what they allow and wanted to implement them in my code. I did that although even at the time of working on it, I could see it wasn't used very well.

## 2024
I was thinking about this project again and the way I implemented the Polynomial class and Matrix class bothered me. Polynomials can be represented in a matrix so having two different objects
which could be more condensed but wasn't and has very specific features to each just felt off. I then had the thought of "Why am I creating classes for these? Numpy probably does it better anyway".
That's when I figured I'd rewrite this code using Numpy and removed those classes. I wanted to create things on my own but it that doesn't always make sense.
Why do that when there are better more complete implementations that can be used. I also cleaned up the code, refactored variable and function names and I have something I feel more proud of.
I added more details about how the function arguments are used and the units of the data. The way I organized my data on the Main.py file feels better, passing in the planets into a list
allows easier control as to what to include in the calculations, if I want to quickly remove a planet or add another and see how things are affected. This is something I'm starting to feel more proud of.
The one thing that does disappoint me is it runs way worse than my original code. That will most likely be the next problem I tackle when I update this project in the future.