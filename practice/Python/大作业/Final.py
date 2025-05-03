## 编写学校成员类，具有成员姓名、年龄和总人数属性，
# 教师类继承学校成员类，具有工号和工资属性，
# 学生类继承学校成员类，具有学号和成绩属性。
# 并实现如下功能：
# （1）创建对象时，总人数加1，对象注销时，总人数减1；
# （2）输出学生的平均成绩。
import gc

class SchoolMember:
    num_total = 0
    num_teacher = 0
    printFlag = True
    def __init__(self, name, age):
        self.name = name
        self.age = age
    def showMember(self):
        if self.printFlag:
            print("当前学校有{}位成员，其中教师{}位，学生{}位".format(self.num_total, self.num_teacher, self.num_total-self.num_teacher))
    @classmethod
    def add_member(cls,isTeacher=False):
        cls.num_total += 1
        if isTeacher:
            cls.num_teacher += 1
    @classmethod
    def del_member(cls,nameIn,isTeacher=False):
        cls.num_total -= 1
        if isTeacher:
            cls.num_teacher -= 1
            if cls.printFlag:
                print("删除教师{}".format(nameIn))
        if not isTeacher:
            if cls.printFlag:
                print("删除学生{}".format(nameIn))

class Teacher(SchoolMember):
    def __init__(self, name, age, id, salary):
        super().__init__(name, age)
        self.id = id
        self.salary = salary
        SchoolMember.add_member(True)
        print("新增教师{}".format(self.name))
        self.showMember()
    def __del__(self):
        SchoolMember.del_member(self.name,True)
        self.showMember()
        pass


class Student(SchoolMember):
    def __init__(self, name, age, id, grade):
        super().__init__(name, age)
        self.id = id
        self.grade = grade
        SchoolMember.add_member()
        print("新增学生{}".format(self.name))
        self.showMember()
    @staticmethod
    def avg_grade(students):
        total_grade = 0
        for student in students:
            total_grade += student.grade
        return total_grade / len(students)
    def __del__(self):
        SchoolMember.del_member(self.name,False)
        self.showMember()
        pass
    def printMeanGrade(self):
        studentList = [obj for obj in gc.get_objects() if isinstance(obj,Student)]
        length = studentList.__len__()
        grade = sum([obj.grade for obj in studentList])
        print("学生平均成绩为{}".format(grade/length))


zhang3 = Teacher("zhang3", 35, "001", 5000)
li = Teacher("li", 40, "002", 6000)
mike = Student("mike", 18, "003", 90)
llh = Teacher("llh", 40, "004", 6000)
licy = Student("licy", 18, "003", 95)
mike.printMeanGrade()
del zhang3
albort = Student("albort", 18, "003", 100)
licy.printMeanGrade()
del li
del licy
albort.printMeanGrade() 
del mike
albort.printMeanGrade()
SchoolMember.printFlag = False